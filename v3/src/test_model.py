"""Independent numerical controls and explicit identifiability checks."""
from pathlib import Path
import sys,time,json
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import norm
from model import *
from prepare import V,dump

def main():
    rng=np.random.default_rng(29103190);n=800;m=8
    X=np.column_stack([np.ones(n),rng.normal(size=n),rng.binomial(1,.5,n)])
    G=rng.binomial(2,.3,(n,m)).astype(float);y=.4*X[:,1]+.15*G[:,0]+rng.normal(size=n)
    d=rng.binomial(1,expit(-1+.2*X[:,1]+.3*G[:,0]))
    c,b,p=marginal(X,y,d,G);err=[]
    for j in range(m):
        A=np.column_stack([X,G[:,j]])
        objective=lambda beta:np.sum(np.logaddexp(0,A@beta)-d*(A@beta))
        fit=minimize(objective,np.zeros(A.shape[1]),jac=lambda beta:A.T@(expit(A@beta)-d),method='BFGS',options={'gtol':1e-8})
        err.append(abs(fit.x[-1]-b['beta'][j]))
        lin=np.linalg.lstsq(A,y,rcond=None)[0]
        assert abs(lin[-1]-c['beta'][j])<1e-12
    assert max(err)<2e-6
    cf,bf,pf=marginal(X,y,d,2-G)
    flip=max(np.max(abs(p[k]-pf[k])) for k in p);assert flip<1e-10
    for j in range(m):assert abs(p['iut'][j]-max(p['p1'][j],p['p2'][j]))<1e-15
    out=coloc(c,b);assert abs(out.sum()-1)<1e-12
    assert np.max(abs(out-coloc(b,c)))>.0 # Trait prior scales differ by design.
    # Exactly equal variant columns cannot identify two different causal configurations.
    h=np.load(V/'reference/region_2.npz')['YRI_study'];meta=json.loads((V/'reference/region_2.json').read_text())['YRI']
    assert np.array_equal(h[:,meta['anchor']],h[:,meta['high']])
    left=generate(np.random.default_rng(123),h,meta,'shared_same','ideal',800)
    right=generate(np.random.default_rng(123),h,meta,'distinct_high','ideal',800)
    assert all(np.array_equal(left[k],right[k]) for k in range(4))
    fixture=V/'results/tests';fixture.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(fixture/'numerical_fixture.npz',X=X,y=y,d=d,G=G,c_beta=c['beta'],c_var=c['model_var'],c_robust=c['var'],b_beta=b['beta'],b_var=b['model_var'],b_robust=b['var'],split_c_beta=c['beta'],split_c_var=c['model_var'],split_b_beta=b['beta'],split_b_var=b['model_var'],pp=out)
    # Time one complete full-panel dataset and the unchanged copula diagnostic.
    hap=np.load(V/'reference/region_0.npz')['YRI_study'];meta=json.loads((V/'reference/region_0.json').read_text())['YRI']
    X,y,d,G,info,_,_=generate(np.random.default_rng(7),hap,meta,'shared_same','realistic')
    t=time.perf_counter();row,_=analyse_dataset(X,y,d,G,meta,True);elapsed=time.perf_counter()-t
    assert row['copula_status']=='ok',row
    result=dict(passed=True,max_logistic_vs_independent_bfgs_error=max(err),max_allele_flip_p_error=flip,perfect_LD_shared_and_distinct_observationally_identical=True,one_dataset_seconds=elapsed,example=row)
    dump(V/'results/tests/numerical_controls.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
