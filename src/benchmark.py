"""Resumable independent-seed calibration and evaluation, with explicit failures."""
from pathlib import Path
import os,json,csv,hashlib,time,argparse,traceback
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import numpy as np
from model import fit,score_scan,full_effects,components,METHODS
from simulate import generate,EFFECTS

ROOT=Path(__file__).resolve().parents[1]
CFG=json.loads((ROOT/'protocol/design.json').read_text())


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def codehash():
    names=['protocol/design.json','src/model.py','src/simulate.py','src/benchmark.py']
    return hashlib.sha256(b''.join((ROOT/n).read_bytes() for n in names)).hexdigest()


def cells():
    out=[];p=CFG['primary'];s=CFG['sensitivity']
    for n in p['n']:
        for env in p['environments']:
            out.append(dict(n=n,environment=env,prevalence=p['prevalence'],maf=p['maf'],ld=0,
                group='primary',setting='primary',cal=p['calibration_null'],null=p['evaluation_null'],
                alt=p['evaluation_per_alternative'],effect_reps=p['full_fit_replicates_per_scenario']))
    for setting in s['settings']:
        for env in s['environments']:
            out.append(dict(n=s['n'],environment=env,prevalence=setting['prevalence'],maf=setting['maf'],ld=setting['ld'],
                group='sensitivity',setting=setting['name'],cal=s['calibration_null'],null=s['evaluation_null'],
                alt=s['evaluation_per_alternative'],effect_reps=s['full_fit_replicates_per_scenario']))
    for i,c in enumerate(out):c['id']=i
    return out


def write_csv(path,rows):
    if not rows:return
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def evaluate(data,with_effects):
    x,w,y,d,g=[data[k] for k in ['x','w','y','d','g']]
    scores={};null_rows=[];effects=[]
    for label,ww in [('Independence',w[:,:0]),('Constant',w[:,:1]),('Varying',w)]:
        started=time.perf_counter()
        try:
            fitted=fit(x,ww,y,d)
            scans=score_scan(fitted,x,ww,y,d,g)
            if label=='Independence':
                targets=[('Marginal Bonferroni','bonf_p'),('Marginal covariance-aware','robust_p'),('Independence likelihood','model_p')]
            else:targets=[(label+' model score','model_p'),(label+' sandwich score','robust_p')]
            for name,key in targets:
                ps=[z.get(key,float('nan')) for z in scans]
                scores[name]=dict(p=ps,status=scans[0]['status'],all_status='ok' if all(z['status']=='ok' for z in scans) else 'panel_failure',
                                  seconds=time.perf_counter()-started)
            rho=components(fitted['theta'],x,ww,y,d)['rho']
            null_rows.append(dict(model=label,status='ok',gradient=fitted['gradient'],
                iterations=fitted['iterations'],starts=fitted['starts'],seconds=time.perf_counter()-started,
                rho_intercept=float(fitted['theta'][2*x.shape[1]+1]) if ww.shape[1] else 0.,
                rho_age=float(fitted['theta'][-1]) if ww.shape[1]==2 else 0.,
                rho_curve_RMSE=float(np.sqrt(np.mean((rho-data['rho'])**2)))))
        except (ValueError,np.linalg.LinAlgError,FloatingPointError) as e:
            names=METHODS[:3] if label=='Independence' else [label+' model score',label+' sandwich score']
            for name in names:scores[name]=dict(p=[float('nan')]*g.shape[1],status=str(e),all_status=str(e),seconds=time.perf_counter()-started)
            null_rows.append(dict(model=label,status=str(e),gradient=float('nan'),iterations=0,starts=0,
                seconds=time.perf_counter()-started,rho_intercept=float('nan'),rho_age=float('nan'),rho_curve_RMSE=float('nan')))
        if with_effects and label!='Independence':
            try:
                ff,betas,ses,model_ses=full_effects(x,ww,y,d,g[:,0])
                for k,trait in enumerate(['continuous','binary_liability']):
                    truth=data['effects'][k]
                    effects.append(dict(model=label,trait=trait,status='ok',truth=truth,estimate=float(betas[k]),
                        sandwich_se=float(ses[k]),model_se=float(model_ses[k]),
                        sandwich_coverage=int(abs(betas[k]-truth)<=1.95996398454*ses[k]),
                        model_coverage=int(abs(betas[k]-truth)<=1.95996398454*model_ses[k]),
                        gradient=ff['gradient']))
            except (ValueError,np.linalg.LinAlgError,FloatingPointError) as e:
                for k,trait in enumerate(['continuous','binary_liability']):
                    effects.append(dict(model=label,trait=trait,status=str(e),truth=data['effects'][k],estimate=float('nan'),
                        sandwich_se=float('nan'),model_se=float('nan'),sandwich_coverage=0,model_coverage=0,gradient=float('nan')))
    return scores,null_rows,effects


def one_cell(cell,smoke=False):
    out=ROOT/('work/smoke' if smoke else 'results/cells')/f"cell_{cell['id']:02d}"
    out.mkdir(parents=True,exist_ok=True)
    identity=codehash()
    if (out/'complete.json').exists() and not smoke:
        old=json.loads((out/'complete.json').read_text())
        assert old['code_sha256']==identity,'Completed-cell source changed; archive old outputs first'
        for name,h in old['outputs'].items():assert digest(out/name)==h
        return dict(cell=cell['id'],reused=True,seconds=0)
    started=time.perf_counter(); calibration=[]
    progress=ROOT/'work/progress';progress.mkdir(parents=True,exist_ok=True)
    def update(phase,scenario,rep,total):
        (progress/f"cell_{cell['id']:02d}.json").write_text(json.dumps(dict(cell=cell['id'],phase=phase,
            scenario=scenario,replicate=rep,total=total,seconds=time.perf_counter()-started,
            utc=datetime.now(timezone.utc).isoformat()))+'\n')
    ncal=3 if smoke else cell['cal']
    for rep in range(ncal):
        dat=generate(99881 if smoke else CFG['seed_calibration'],cell['id'],rep,'null',cell,CFG['markers'])
        ss,_,_=evaluate(dat,False)
        for method,z in ss.items():calibration.append(dict(replicate=rep,method=method,p=z['p'][0],status=z['status']))
        if rep%25==0:update('calibration','null',rep,ncal)
    write_csv(out/'calibration.csv',calibration)
    thresholds={};cal_fail={}
    for method in METHODS:
        rows=[r for r in calibration if r['method']==method]
        cal_fail[method]=sum(r['status']!='ok' or not np.isfinite(r['p']) for r in rows)
        ps=np.sort([r['p'] for r in rows])
        rank=max(1,int(np.floor(.05*(ncal+1))))
        thresholds[method]=float(ps[rank-1]) if not cal_fail[method] else None
    freeze=dict(code_sha256=identity,calibration_sha256=digest(out/'calibration.csv'),thresholds=thresholds,
                calibration_failures=cal_fail,frozen_before_evaluation_utc=datetime.now(timezone.utc).isoformat())
    (out/'threshold_freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
    score_rows=[];fit_rows=[];effect_rows=[];dataset_rows=[]
    for scenario in EFFECTS:
        nrep=2 if smoke else cell['null'] if scenario=='null' else cell['alt']
        for rep in range(nrep):
            dat=generate(99882 if smoke else CFG['seed_evaluation'],cell['id'],rep,scenario,cell,CFG['markers'])
            ss,ff,ee=evaluate(dat,rep<(1 if smoke else cell['effect_reps']))
            uid=hashlib.sha256(dat['g'].tobytes()+dat['y'].tobytes()+dat['d'].tobytes()).hexdigest()
            common=dict(cell=cell['id'],scenario=scenario,replicate=rep)
            dataset_rows.append(dict(**common,fingerprint=uid,cases=int(dat['d'].sum()),n=cell['n'],
                binary_intercept=dat['binary_intercept']))
            for method,z in ss.items():
                ps=z['p'];p=ps[0];cut=thresholds[method]
                score_rows.append(dict(**common,method=method,status=z['status'],panel_status=z['all_status'],
                    p=p,nominal_05=int(p<.05),nominal_01=int(p<.01),
                    calibrated_05=int(p<cut) if cut is not None else '',
                    panel_FWER_05=int(any(q<.05/CFG['markers'] for q in ps)),
                    threshold=cut if cut is not None else '',seconds=z['seconds'],
                    **{'p'+str(j):q for j,q in enumerate(ps)}))
            fit_rows.extend(dict(**common,**row) for row in ff)
            effect_rows.extend(dict(**common,**row) for row in ee)
            if rep%25==0:update('evaluation',scenario,rep,nrep)
    write_csv(out/'scores.csv',score_rows);write_csv(out/'null_fits.csv',fit_rows)
    write_csv(out/'effects.csv',effect_rows);write_csv(out/'datasets.csv',dataset_rows)
    outputs={p.name:digest(p) for p in out.iterdir() if p.suffix in ['.csv','.json'] and p.name!='complete.json'}
    done=dict(cell=cell,code_sha256=identity,outputs=outputs,
        evaluation_datasets=len(dataset_rows),score_records=len(score_rows),
        seconds=time.perf_counter()-started,completed_utc=datetime.now(timezone.utc).isoformat())
    (out/'complete.json').write_text(json.dumps(done,indent=2)+'\n')
    update('complete','all',len(dataset_rows),len(dataset_rows))
    return dict(cell=cell['id'],seconds=done['seconds'],datasets=len(dataset_rows))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--smoke',action='store_true');args=ap.parse_args()
    controls=json.loads((ROOT/'results/numerical_controls.json').read_text())
    cross=json.loads((ROOT/'results/R_validation.json').read_text())
    assert controls['passed'] and controls['engine_sha256']==digest(ROOT/'src/model.py')
    assert cross['passed'] and cross['engine_md5']==hashlib.md5((ROOT/'src/model.py').read_bytes()).hexdigest()
    planned=cells()
    (ROOT/'protocol/cells.json').write_text(json.dumps(planned,indent=2)+'\n')
    if args.smoke:
        for c in [planned[2],planned[3]]:print(one_cell(c,True),flush=True)
        return
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        pending={pool.submit(one_cell,c):c for c in planned}
        for f in as_completed(pending):print(json.dumps(f.result()),flush=True)

if __name__=='__main__':main()
