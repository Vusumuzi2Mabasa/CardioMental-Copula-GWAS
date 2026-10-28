"""Chunked benchmark with calibration barrier, frozen sources and explicit failures."""
from pathlib import Path
import argparse,csv,json,hashlib,time,os,sys
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import numpy as np
from oracle import METHODS,evaluate
from generate import generate

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parent
CORE=['src/model.py','src/simulate.py','v2/src/oracle.py','v2/src/generate.py','v2/src/run.py','v2/protocol/design.json']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def identity():return hashlib.sha256(b''.join((REPO/p).read_bytes() for p in CORE)).hexdigest()
def config():return json.loads((ROOT/'protocol/design.json').read_text())
def cells():
    return [dict(id=i,n=n,environment=e,maf=.25,prevalence=.15) for i,(n,e) in enumerate((n,e) for n in [800,1600] for e in ['constant','varying'])]
def scenarios(cell):
    out=[dict(name='null',beta=[0.,0.],nrep=config()['evaluation_null'])]
    for level,scale in zip(['weak','moderate'],config()['scales']):
        for arch,b in config()['architectures'].items():
            beta=(np.array(b)*scale*np.sqrt(800/cell['n'])).tolist()
            out.append(dict(name=level+'_'+arch,beta=beta,nrep=config()['evaluation_per_alternative']))
    return out
def fingerprint(dat):return hashlib.sha256(dat['g'].tobytes()+dat['y'].tobytes()+dat['d'].tobytes()).hexdigest()
def read(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))
def write(path,rows):
    with Path(path).open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
def chunk(cell,phase,scenario,sid,start,stop,smoke):
    base=ROOT/('work/smoke' if smoke else 'results/chunks')/f"cell_{cell['id']}"
    base.mkdir(parents=True,exist_ok=True)
    name=f'{phase}_{sid}_{start:05d}';out=base/(name+'.csv');done=base/(name+'.json')
    code=identity()
    if done.exists():
        old=json.loads(done.read_text());assert old['code_sha256']==code and old['csv_sha256']==sha(out)
        return old
    began=time.perf_counter();rows=[];cfg=config()
    seed=cfg['seed_calibration'] if phase=='calibration' else cfg['seed_evaluation']
    if smoke:seed+=555555
    for rep in range(start,stop):
        dat=generate(seed,cell,rep,scenario['beta'],stream=sid)
        ss=evaluate(dat);uid=fingerprint(dat)
        for method in METHODS:
            z=ss[method];rows.append(dict(cell=cell['id'],phase=phase,scenario=scenario['name'],replicate=rep,
                method=method,fingerprint=uid,cases=int(dat['d'].sum()),beta1=scenario['beta'][0],beta2=scenario['beta'][1],**z))
    write(out,rows)
    result=dict(cell=cell['id'],phase=phase,scenario=scenario['name'],start=start,stop=stop,rows=len(rows),
                failures=sum(r['status']!='ok' for r in rows),code_sha256=code,csv_sha256=sha(out),
                seconds=time.perf_counter()-began,completed_utc=datetime.now(timezone.utc).isoformat())
    done.write_text(json.dumps(result,indent=2)+'\n')
    return result

def freeze(smoke=False):
    cfg=config();base=ROOT/('work/smoke' if smoke else 'results/chunks');frozen=[]
    ncal=4 if smoke else cfg['calibration_null']
    for cell in cells():
        rows=[r for p in sorted((base/f"cell_{cell['id']}").glob('calibration_*.csv')) for r in read(p)]
        assert len(rows)==ncal*len(METHODS)
        for method in METHODS:
            rr=[r for r in rows if r['method']==method];assert sorted(int(r['replicate']) for r in rr)==list(range(ncal))
            failures=sum(r['status']!='ok' or not np.isfinite(float(r['p'])) for r in rr)
            for alpha in cfg['alpha_grid']:
                rank=max(1,int(np.floor(alpha*(ncal+1))))
                cutoff=float(np.sort([float(r['p']) for r in rr])[rank-1]) if failures==0 else None
                frozen.append(dict(cell=cell['id'],method=method,alpha=alpha,rank=rank,ncal=ncal,
                    target=rank/(ncal+1),cutoff=cutoff,failures=failures))
    out=ROOT/('work/smoke/thresholds.json' if smoke else 'results/thresholds.json')
    value=dict(code_sha256=identity(),frozen_before_evaluation_utc=datetime.now(timezone.utc).isoformat(),thresholds=frozen,
               calibration_hashes={str(p.relative_to(ROOT)):sha(p) for p in sorted(base.glob('cell_*/calibration_*.csv'))})
    if out.exists():
        old=json.loads(out.read_text());assert old['thresholds']==frozen and old['calibration_hashes']==value['calibration_hashes']
        return
    out.write_text(json.dumps(value,indent=2)+'\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=2);ap.add_argument('--smoke',action='store_true');args=ap.parse_args()
    assert json.loads((ROOT/'results/numerical_controls.json').read_text())['passed']
    assert json.loads((ROOT/'results/R_validation.json').read_text())['passed']
    if not args.smoke:
        lock=json.loads((ROOT/'protocol/source_lock.json').read_text());assert lock['code_sha256']==identity()
        for path,h in lock['files'].items():assert sha(REPO/path)==h,path
    cfg=config();size=cfg['chunk_size'];began=time.perf_counter()
    for phase in ['calibration','evaluation']:
        tasks=[]
        for cell in cells():
            ss=[dict(name='null',beta=[0.,0.],nrep=4 if args.smoke else cfg['calibration_null'])] if phase=='calibration' else scenarios(cell)
            for sid,sc in enumerate(ss):
                nrep=2 if args.smoke and phase=='evaluation' else sc['nrep']
                for first in range(0,nrep,size):tasks.append((cell,phase,sc,sid,first,min(first+size,nrep),args.smoke))
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            ff=[pool.submit(chunk,*t) for t in tasks]
            for i,f in enumerate(as_completed(ff)):
                z=f.result()
                if i%10==0 or i==len(ff)-1:print(json.dumps(dict(phase=phase,completed_chunks=i+1,total_chunks=len(ff),elapsed_seconds=time.perf_counter()-began,last=z)),flush=True)
        if phase=='calibration':freeze(args.smoke)
    print('Benchmark complete',time.perf_counter()-began,flush=True)
if __name__=='__main__':main()
