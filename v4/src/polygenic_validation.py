"""Uniform larger-bank calibration and fresh cohorts, preserving the initial experiment."""
from pathlib import Path
import sys,json,csv,time,os,argparse,hashlib,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.append(str(ROOT/'v3/src'))
from prepare import dump,table,sha
from polygenic import simulate,tag,V
from kernel import kernel
from reml import fit
from freeze import verify
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--freeze',action='store_true');args=ap.parse_args()
    original=verify();spec=V/'protocol/validation_extension.json';me=Path(__file__)
    identity=hashlib.sha256((original['identity']+sha(spec)+sha(me)).encode()).hexdigest()
    lock=V/'protocol/validation_source_lock.json'
    if args.freeze:
        assert not lock.exists()
        dump(lock,dict(identity=identity,original_identity=original['identity'],driver_sha256=sha(me),extension_sha256=sha(spec),parent_commit=subprocess.check_output(['git','-c','core.fsmonitor=false','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),stage='Before all fresh validation datasets'))
        print(identity);return
    assert json.loads(lock.read_text())['identity']==identity
    base=V/'results/polygenic_validation';out=base/'polygenic_chunks';out.mkdir(parents=True,exist_ok=True)
    for phase,count,offset in [('cal',1000,1000000),('eval',200,2000000)]:
        if phase=='eval':assert (base/'polygenic_thresholds.json').exists()
        for pop in ['YRI','LWK']:
            for n in [600,1200]:
                with np.load(kernel(pop,n)) as stored:panel={k:stored[k] for k in stored.files}
                lam=panel['eigenvalues']
                for re in [0.,.4]:
                    cell=tag(pop,n,re)
                    for rg in ([0.] if phase=='cal' else [-.5,0.,.5]):
                        for start in range(0,count,10):
                            dest=out/f'{phase}_{cell}_g{rg:.1f}_{start:04d}.csv';meta=dest.with_suffix('.json')
                            if dest.exists() and meta.exists():
                                m=json.loads(meta.read_text());assert m['identity']==identity and m['sha256']==sha(dest);continue
                            records=[];began=time.perf_counter()
                            for rep in range(start,min(start+10,count)):
                                latent,observed,info=simulate(panel,pop,n,re,rg,phase,offset+rep)
                                for mode,y in [('liability',latent),('observed_binary',observed)]:
                                    row=dict(population=pop,n=n,residual_truth=re,genetic_truth=rg,phase=phase,replicate=rep,seed_replicate=offset+rep,mode=mode,**info)
                                    try:
                                        stats,theta,tn=fit(lam,y);row.update(stats)
                                        if phase=='eval' and rep==0:
                                            fixture=base/'polygenic_fixtures';fixture.mkdir(parents=True,exist_ok=True)
                                            dump(fixture/f'{cell}_g{rg:.1f}_{mode}.json',dict(eigenvalues=lam.tolist(),y=y.tolist(),theta=theta.tolist(),theta_null=tn.tolist(),expected=stats))
                                    except (ValueError,np.linalg.LinAlgError) as e:
                                        row.update(status=str(e),genetic_c=0.,genetic_m=0.,genetic_covariance=0.,genetic_correlation=0.,residual_c=0.,residual_m=0.,residual_covariance=0.,residual_correlation=0.,LR=0.,asymptotic_p=1.,lower=-1.,upper=1.,ci_status='failed',boundary=0,max_gradient=0.,loglik_alt=0.,loglik_null=0.)
                                    records.append(row)
                            tmp=dest.with_suffix('.tmp');table(tmp,records);os.replace(tmp,dest)
                            dump(meta,dict(identity=identity,sha256=sha(dest),rows=len(records),seconds=time.perf_counter()-began))
                        print('fresh',phase,cell,'rg',rg,'complete',flush=True)
        if phase=='cal':freeze(base)
def freeze(base):
    banks={}
    for p in (base/'polygenic_chunks').glob('cal_*.csv'):
        with p.open() as f:
            for r in csv.DictReader(f):banks.setdefault(tag(r['population'],int(r['n']),float(r['residual_truth']))+'_'+r['mode'],[]).append(r)
    assert len(banks)==16
    values={}
    for cell,rows in banks.items():
        assert len(rows)==1000;fail=sum(r['status']!='ok' for r in rows)
        values[cell]=dict(active=not fail,failures=fail,LR_cutoff=float(np.sort([float(r['LR']) for r in rows])[950]))
    path=base/'polygenic_thresholds.json';obj=dict(thresholds=values,strict=True,alpha=.05,calibration_datasets=8000,outcome_modes_not_independent=True)
    if path.exists():assert json.loads(path.read_text())==obj
    else:dump(path,obj)
    print('Fresh validation thresholds frozen',flush=True)
if __name__=='__main__':main()
