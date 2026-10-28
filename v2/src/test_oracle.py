"""Independent expanded-likelihood finite differences and oracle controls."""
from pathlib import Path
import json,hashlib,time
import numpy as np
from scipy.optimize._numdiff import approx_derivative
from scipy.stats import chi2
from oracle import components,restricted_fit,restricted_score,true_theta,evaluate
from generate import generate

ROOT=Path(__file__).resolve().parents[1]
def run():
    began=time.perf_counter();cell=dict(id=0,n=800,maf=.25,prevalence=.15,environment='varying')
    data=generate(487119,cell,0,(0.,0.));x,w,y,d,g=[data[k] for k in ['x','w','y','d','g']]
    th=restricted_fit(data);err=[];references=[]
    for allknown in [False,True]:
        theta=true_theta(data) if allknown else th
        gg=g[:,0]-(.5 if allknown else g[:,0].mean())
        xx=np.column_stack([x,gg]);enlarged=np.r_[theta[:3],0.,theta[3:6],0.,theta[6:]]
        bidx=[3,7]
        # Explicit free indices are SNP effects then seven marginal nuisances.
        active=bidx if allknown else bidx+[0,1,2,4,5,6,8]
        def scores(v):
            tt=enlarged.copy();tt[active]=v
            return components(tt,xx,w,y,d)['scores'][:,active]
        v=enlarged[active];ss=scores(v)
        h=-approx_derivative(lambda a:scores(a).sum(0),v)
        u=ss[:,:2].sum(0);info=h[:2,:2];psi=ss[:,:2]
        if not allknown:
            a=h[:2,2:]@np.linalg.inv(h[2:,2:]);u-=a@ss[:,2:].sum(0)
            info=info-a@h[2:,:2];psi=psi-ss[:,2:]@a.T
        tm=float(u@np.linalg.solve(info,u));tr=float(u@np.linalg.solve(psi.T@psi,u))
        direct=restricted_score(theta,data,allknown)
        err.extend([abs(tm-direct['model_stat']),abs(tr-direct['robust_stat'])])
        assert max(err)<1e-4,err
        # Allele coding invariance includes inversion of the known population mean.
        flipped=dict(data,g=2-g,maf=1-data['maf'])
        z=restricted_score(theta,flipped,allknown)
        assert abs(z['model_p']-direct['model_p'])<1e-10
        references.append(dict(all_known=allknown,theta=theta.tolist(),**direct))
    baseline=evaluate(data);assert all(v['status']=='ok' for v in baseline.values())
    # The data generator must not recalibrate the intercept under alternatives.
    alternative=generate(487119,cell,0,(.05,.08))
    assert alternative['binary_intercept']==data['binary_intercept']
    independent=dict(data,alpha=np.zeros(2))
    theta=true_theta(independent)
    cc=components(theta,x,w,y,d)
    from scipy.stats import norm
    wanted=norm.logpdf(y,loc=x@theta[:3])+norm.logcdf((2*d-1)*(x@theta[3:6]))
    assert np.max(abs(wanted-cc['loglik']))<1e-10
    np.savetxt(ROOT/'results/oracle_fixture.csv',np.column_stack([x[:,1:],y,d,g]),delimiter=',',
               header='age,exposure,y,d,g',comments='')
    (ROOT/'results/oracle_fixture_expected.json').write_text(json.dumps(references,indent=2)+'\n')
    result=dict(passed=True,expanded_hessian_max_stat_error=max(err),allele_flip_passed=True,
                fixed_intercept_passed=True,independence_passed=True,seconds=time.perf_counter()-began)
    (ROOT/'results/numerical_controls.json').write_text(json.dumps(result,indent=2)+'\n');print(result)

if __name__=='__main__':run()
