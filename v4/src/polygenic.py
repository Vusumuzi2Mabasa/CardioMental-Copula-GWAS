"""Independent calibration and evaluation of genome-scale synthetic genetic correlation."""
from pathlib import Path
import sys,json,argparse,csv,time,os,hashlib
import numpy as np
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[2];sys.path.append(str(ROOT/'v3/src'))
from prepare import dump,table,sha
from kernel import kernel,V
from reml import fit

def simulate(panel,pop,n,re,rg,phase,rep):
    seed=[29103111 if phase=='cal' else 29103112,['YRI','LWK'].index(pop),n,int(re*10),int((rg+1)*10),rep]
    rng=np.random.default_rng(np.random.SeedSequence(seed));Z,B,lam=[panel[k] for k in ['Z','B','eigenvalues']]
    gchol=np.linalg.cholesky(.4*np.array([[1,rg],[rg,1]]));echol=np.linalg.cholesky(.6*np.array([[1,re],[re,1]]))
    effects=rng.normal(size=(Z.shape[1],2))@gchol.T
    y=Z@effects/np.sqrt(Z.shape[1])+rng.normal(size=(n,2))@echol.T
    modes=B.T@y;truth=(y[:,1]>norm.isf(.15)).astype(int)
    d=np.where(truth,rng.uniform(size=n)<.9,rng.uniform(size=n)>.98).astype(float)
    observed=np.column_stack([modes[:,0],B.T@d])
    return modes,observed,dict(prevalence=float(d.mean()),phenotypic_correlation=float(np.corrcoef(y[:,0],d)[0,1]),liability_phenotypic_correlation=float(np.corrcoef(y.T)[0,1]))
def tag(pop,n,re):return f'{pop}_n{n}_e{re:.1f}'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['cal','eval','all'],default='all');ap.add_argument('--smoke',action='store_true');args=ap.parse_args()
    out=V/('work/polygenic_smoke' if args.smoke else 'results/polygenic_chunks');out.mkdir(parents=True,exist_ok=True)
    from freeze import verify
    lock=verify()['identity']
    for phase in (['cal','eval'] if args.phase=='all' else [args.phase]):
        if phase=='eval' and not args.smoke:
            assert (V/'results/polygenic_thresholds.json').exists(),'Complete and freeze calibration first'
        for pop in (['YRI'] if args.smoke else ['YRI','LWK']):
            for n in ([600] if args.smoke else [600,1200]):
                with np.load(kernel(pop,n)) as stored:panel={k:stored[k] for k in stored.files}
                lam=panel['eigenvalues']
                for re in ([0.] if args.smoke else [0.,.4]):
                    cell=tag(pop,n,re)
                    for rg in ([0.] if phase=='cal' else [-.5,0.,.5]):
                        count=1 if args.smoke else (200 if phase=='cal' else 100)
                        for start in range(0,count,10):
                            dest=out/f'{phase}_{cell}_g{rg:.1f}_{start:03d}.csv';meta=dest.with_suffix('.json')
                            if dest.exists() and meta.exists():
                                d=json.loads(meta.read_text());assert d['identity']==lock and d['sha256']==sha(dest);continue
                            records=[];t=time.perf_counter()
                            for rep in range(start,min(start+10,count)):
                                latent,observed,info=simulate(panel,pop,n,re,rg,phase,rep)
                                for mode,y in [('liability',latent),('observed_binary',observed)]:
                                    row=dict(population=pop,n=n,residual_truth=re,genetic_truth=rg,phase=phase,replicate=rep,mode=mode,**info)
                                    try:
                                        stats,theta,tn=fit(lam,y);row.update(stats)
                                        if phase=='eval' and rep==0:
                                            fixture=V/('work/polygenic_smoke_fixtures' if args.smoke else 'results/polygenic_fixtures');fixture.mkdir(parents=True,exist_ok=True)
                                            dump(fixture/f'{cell}_g{rg:.1f}_{mode}.json',dict(eigenvalues=lam.tolist(),y=y.tolist(),theta=theta.tolist(),theta_null=tn.tolist(),expected=stats))
                                    except (ValueError,np.linalg.LinAlgError) as e:
                                        row.update(status=str(e),genetic_c=0.,genetic_m=0.,genetic_covariance=0.,genetic_correlation=0.,residual_c=0.,residual_m=0.,residual_covariance=0.,residual_correlation=0.,LR=0.,asymptotic_p=1.,lower=-1.,upper=1.,ci_status='failed',boundary=0,max_gradient=0.,loglik_alt=0.,loglik_null=0.)
                                    records.append(row)
                            tmp=dest.with_suffix('.tmp');table(tmp,records);os.replace(tmp,dest)
                            dump(meta,dict(identity=lock,sha256=sha(dest),rows=len(records),seconds=time.perf_counter()-t))
                        print(phase,cell,'rg',rg,'complete',flush=True)
        if phase=='cal' and not args.smoke:freeze_thresholds()
def freeze_thresholds():
    banks={}
    for p in (V/'results/polygenic_chunks').glob('cal_*.csv'):
        with p.open() as f:
            for r in csv.DictReader(f):banks.setdefault(tag(r['population'],int(r['n']),float(r['residual_truth']))+'_'+r['mode'],[]).append(r)
    assert len(banks)==16
    values={}
    for cell,rows in banks.items():
        assert len(rows)==200
        fail=sum(r['status']!='ok' for r in rows)
        values[cell]=dict(active=not fail,failures=fail,LR_cutoff=float(np.sort([float(r['LR']) for r in rows])[190]))
    dest=V/'results/polygenic_thresholds.json';obj=dict(thresholds=values,strict=True,alpha=.05,calibration_datasets=1600,outcome_modes_not_independent=True)
    if dest.exists():assert json.loads(dest.read_text())==obj
    else:dump(dest,obj)
    print('Polygenic thresholds frozen',flush=True)
if __name__=='__main__':main()
