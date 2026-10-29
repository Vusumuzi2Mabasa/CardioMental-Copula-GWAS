"""Frozen, resumable simulation with independent calibration and evaluation seeds."""
from pathlib import Path
import argparse, csv, hashlib, json, os, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from model import generate,analyse_dataset,SCENARIOS,NO_SHARED
from prepare import V,ROOT,sha,dump,table

def identity():
    paths=['v3/protocol/design.json','v3/src/model.py','v3/src/run.py']
    paths += [str(p.relative_to(ROOT)).replace('\\','/') for p in sorted((V/'reference').glob('*')) if p.is_file()]
    return hashlib.sha256(json.dumps({p:sha(ROOT/p) for p in paths},sort_keys=True).encode()).hexdigest()
def cellid(pop,region,regime):return f'{pop}_r{region}_{regime}'
def task(job):
    pop,region,regime,phase,scenario,start,count,smoke,fingerprint=job
    tag=cellid(pop,region,regime);out=V/('work/smoke' if smoke else 'results/chunks')/tag
    out.mkdir(parents=True,exist_ok=True);dest=out/f'{phase}_{scenario}_{start:04d}.csv';meta=dest.with_suffix('.json')
    if dest.exists() and meta.exists():
        state=json.loads(meta.read_text())
        if state['identity']!=fingerprint or sha(dest)!=state['sha256']:raise ValueError('Checkpoint identity mismatch')
        return 0
    hap=np.load(V/f'reference/region_{region}.npz')[pop+'_study'];causal=json.loads((V/f'reference/region_{region}.json').read_text())[pop]
    records=[];t0=time.perf_counter()
    for rep in range(start,start+count):
        seed=[29103101 if phase=='cal' else 29103102,['YRI','LWK'].index(pop),region,['ideal','realistic'].index(regime),SCENARIOS.index(scenario),rep]
        rng=np.random.default_rng(np.random.SeedSequence(seed))
        X,y,d,G,info,bc,bm=generate(rng,hap,causal,scenario,regime)
        unique=hashlib.sha256(X.tobytes()+y.tobytes()+d.tobytes()+G.tobytes()).hexdigest()
        row=dict(population=pop,region=region,regime=regime,phase=phase,scenario=scenario,replicate=rep,fingerprint=unique,**info)
        try:
            # Fixed first 20 evaluation replicates per cell/architecture form the copula diagnostic panel.
            stats,fixture=analyse_dataset(X,y,d,G,causal,copula=(phase=='eval' and rep<20))
            row.update(status='ok',**stats)
            if phase=='eval' and rep==0 and scenario in ['null','shared_same','distinct_high']:
                fixtures=V/('work/smoke_fixtures' if smoke else 'results/fixtures');fixtures.mkdir(parents=True,exist_ok=True)
                np.savez_compressed(fixtures/f'{tag}_{scenario}.npz',**fixture)
        except (ValueError,np.linalg.LinAlgError,FloatingPointError) as e:
            row.update(status=str(e),iut_p=1.,joint_p=1.,anchor_iut=1.,anchor_joint=1.,split_pp4=0.,split_pp3=0.,overlap_pp4=0.,prior_low_pp4=0.,prior_high_pp4=0.,min_p_cardio=1.,min_p_mental=1.,copula_status='dataset_failed',copula_p=1.)
        records.append(row)
    tmp=dest.with_suffix('.tmp');table(tmp,records);os.replace(tmp,dest)
    dump(meta,dict(identity=fingerprint,sha256=sha(dest),rows=len(records),seconds=time.perf_counter()-t0))
    return len(records)

def thresholds(fingerprint):
    cfg=json.loads((V/'protocol/design.json').read_text());out={}
    for pop in cfg['populations']:
        for region in range(3):
            for regime in cfg['regimes']:
                tag=cellid(pop,region,regime);banks={}
                for scenario in NO_SHARED:
                    rows=[]
                    for path in sorted((V/'results/chunks'/tag).glob(f'cal_{scenario}_*.csv')):
                        with path.open() as f:rows.extend(csv.DictReader(f))
                    if len(rows)!=cfg['calibration_replicates']:raise ValueError('Incomplete calibration')
                    if any(r['status']!='ok' for r in rows):raise ValueError('Calibration fit failure blocks thresholds')
                    banks[scenario]=rows
                def cutoff(key,scenarios,upper):
                    # Strict rejection, using exchangeable order-statistic rank.
                    n=cfg['calibration_replicates'];k=int(np.floor(.05*(n+1)))
                    values=[np.sort([float(r[key]) for r in banks[s]])[n-k if upper else k-1] for s in scenarios]
                    return float(max(values) if upper else min(values))
                out[tag]=dict(split_pp4=cutoff('split_pp4',NO_SHARED,True),overlap_pp4=cutoff('overlap_pp4',NO_SHARED,True),iut_p=cutoff('iut_p',NO_SHARED[:3],False),joint_p=cutoff('joint_p',['null'],False))
    dump(V/'results/thresholds.json',dict(identity=fingerprint,alpha=.05,strict=True,composite_null=NO_SHARED,thresholds=out))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=2);ap.add_argument('--smoke',action='store_true');ap.add_argument('--phase',choices=['cal','eval','all'],default='all');args=ap.parse_args()
    cfg=json.loads((V/'protocol/design.json').read_text());fp=identity();lock=json.loads((V/'protocol/source_lock.json').read_text())
    if fp!=lock['identity']:raise ValueError('Frozen source/reference identity mismatch')
    for phase in (['cal','eval'] if args.phase=='all' else [args.phase]):
        if phase=='eval' and not args.smoke:
            saved=json.loads((V/'results/thresholds.json').read_text())
            if saved['identity']!=fp:raise ValueError('Threshold identity mismatch')
        scenarios=NO_SHARED if phase=='cal' else SCENARIOS
        jobs=[];n=1 if args.smoke else cfg['calibration_replicates' if phase=='cal' else 'evaluation_replicates'];chunk=1 if args.smoke else 10
        for pop in (['YRI'] if args.smoke else cfg['populations']):
            for region in ([0] if args.smoke else range(3)):
                for regime in cfg['regimes']:
                    for scenario in scenarios:
                        for start in range(0,n,chunk):jobs.append((pop,region,regime,phase,scenario,start,min(chunk,n-start),args.smoke,fp))
        completed=0;t=time.perf_counter()
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for done in as_completed([pool.submit(task,j) for j in jobs]):
                completed+=1;new=done.result()
                if completed%20==0 or completed==len(jobs):print(f'{phase}: chunks {completed}/{len(jobs)}, seconds {time.perf_counter()-t:.1f}',flush=True)
        if phase=='cal' and not args.smoke:thresholds(fp)
if __name__=='__main__':main()
