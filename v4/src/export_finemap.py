"""Recreate prespecified v3 subsamples for ancestry-aware multi-signal fine-mapping."""
from pathlib import Path
import sys,json,argparse
import numpy as np
from scipy.special import expit
ROOT=Path(__file__).resolve().parents[2];sys.path.append(str(ROOT/'v3/src'))
from model import generate,linear,logistic,SCENARIOS,NO_SHARED
from prepare import dump,sha
V=ROOT/'v4'
def residual_ld(X,G,weights=None):
    w=np.ones(len(G)) if weights is None else weights
    gr=G-X@np.linalg.solve(X.T@(w[:,None]*X),X.T@(w[:,None]*G))
    cov=gr.T@(w[:,None]*gr);s=np.sqrt(np.diag(cov));R=cov/s[:,None]/s[None,:]
    return .95*R+.05*np.eye(G.shape[1])
def null_weights(X,d):
    b=np.zeros(X.shape[1]);b[0]=np.log(d.mean()/(1-d.mean()))
    for _ in range(50):
        mu=expit(X@b);w=mu*(1-mu);step=np.linalg.solve(X.T@(w[:,None]*X),X.T@(d-mu));b+=step
        if max(abs(step))<1e-10:break
    return expit(X@b)*(1-expit(X@b))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['cal','eval'],required=True);ap.add_argument('--smoke',action='store_true');args=ap.parse_args()
    from freeze import verify
    verify()
    cfg=json.loads((V/'protocol/design.json').read_text())['level2'];out=V/('work/fm_smoke' if args.smoke else 'work/finemap_inputs');out.mkdir(parents=True,exist_ok=True)
    if args.phase=='eval' and not (ROOT/'v3/results/thresholds.json').exists():raise ValueError('Wait for the original v3 calibration freeze')
    if args.phase=='eval' and not args.smoke and not (V/'results/finemap_thresholds.json').exists():raise ValueError('Fine-mapping thresholds must be frozen first')
    n=1 if args.smoke else 20
    for pi,pop in enumerate(['YRI'] if args.smoke else cfg['populations']):
        for region in ([0] if args.smoke else cfg['regions']):
            panel=np.load(ROOT/f'v3/reference/region_{region}.npz');h=panel[pop+'_study'];ref=panel[pop+'_reference'];gref=ref[::2]+ref[1::2]
            Rref=residual_ld(np.ones((len(gref),1)),gref)
            meta=json.loads((ROOT/f'v3/reference/region_{region}.json').read_text())[pop]
            for ri,regime in enumerate(cfg['regimes']):
                for sc in (NO_SHARED if args.phase=='cal' else SCENARIOS):
                    path=out/f'{args.phase}_{pop}_r{region}_{regime}_{sc}.json'
                    if path.exists():continue
                    records=[]
                    for rep in range(n):
                        seed=[29103101 if args.phase=='cal' else 29103102,['YRI','LWK'].index(pop),region,ri,SCENARIOS.index(sc),rep]
                        X,y,d,G,info,bc,bm=generate(np.random.default_rng(np.random.SeedSequence(seed)),h,meta,sc,regime)
                        mid=len(y)//2;xc,yc,gc=X[:mid],y[:mid],G[:mid];xb,db,gb=X[mid:],d[mid:],G[mid:]
                        c=linear(xc,yc,gc);b=logistic(xb,db,gb)
                        records.append(dict(replicate=rep,Nc=len(yc),Nb=len(db),c_beta=c['beta'].tolist(),c_var=c['model_var'].tolist(),b_beta=b['beta'].tolist(),b_var=b['model_var'].tolist(),Rc=residual_ld(xc,gc).tolist(),Rb=residual_ld(xb,gb,null_weights(xb,db)).tolist(),true_c=np.flatnonzero(bc).tolist(),true_b=np.flatnonzero(bm).tolist()))
                    dump(path,dict(population=pop,region=region,regime=regime,scenario=sc,phase=args.phase,Rref=Rref.tolist(),records=records))
                    print(path.name,flush=True)
if __name__=='__main__':main()
