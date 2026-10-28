"""Independent finite differences, block information and null limiting controls."""
from pathlib import Path
import json,time,hashlib
import numpy as np
from scipy.optimize._numdiff import approx_derivative
from scipy.integrate import quad
from scipy.special import ndtr
from scipy.stats import norm,chi2
from model import components,objective,fit,derivatives,score_scan,full_effects
from simulate import generate

ROOT=Path(__file__).resolve().parents[1]


def run():
    began=time.perf_counter()
    cell=dict(n=650,maf=.25,prevalence=.15,environment='varying')
    dat=generate(987654,0,0,'shared',cell)
    x,w,y,d,g=[dat[k] for k in ['x','w','y','d','g']]
    theta=np.array([.1,.25,-.1,-1.,.2,.3,.1,.3,.45])
    f,grad=objective(theta,x,w,y,d)
    numeric=approx_derivative(lambda th:objective(th,x,w,y,d)[0],theta,method='3-point').ravel()
    gradient_error=float(np.max(np.abs(grad-numeric)))
    assert gradient_error<2e-7,gradient_error
    for rho in [-.8,0,.7]:
        for eta in [-2.,.3,1.5]:
            mass=sum(quad(lambda z: norm.pdf(z)*ndtr((2*dd-1)*(eta+rho*z)/np.sqrt(1-rho*rho)),
                          -10,10,epsabs=1e-11)[0] for dd in [0,1])
            assert abs(mass-1)<1e-9
    independence=np.r_[theta[:7],0.,0.]
    ci=components(independence,x,w,y,d)
    expected=norm.logpdf(y,loc=x@theta[:3],scale=np.exp(theta[6]))+norm.logcdf((2*d-1)*(x@theta[3:6]))
    independent_error=float(np.max(np.abs(ci['loglik']-expected)))
    assert independent_error<1e-10
    fits=[];fixture_checks=[]
    score_differences=[]; flip_errors=[]; h_errors=[]
    for ww in [w[:,:0],w[:,:1],w]:
        best=fit(x,ww,y,d,multistart=True)
        ordinary=fit(x,ww,y,d)
        assert abs(best['nll']-ordinary['nll'])<1e-4
        p0=score_scan(best,x,ww,y,d,g)
        pf=score_scan(best,x,ww,y,d,2-g)
        assert all(z['status']=='ok' for z in p0+pf)
        flip_errors.extend(abs(a['robust_p']-b['robust_p']) for a,b in zip(p0,pf))
        h,_,asym=derivatives(best['theta'],x,ww,y,d)
        hn=approx_derivative(lambda th:objective(th,x,ww,y,d)[1]*len(y),best['theta'])
        h_errors.append(float(np.max(np.abs(h-hn))))
        assert np.max(np.abs(h-hn))<2e-4
        # Joint observed Hessian with beta=0 must reproduce Schur complement score.
        p=x.shape[1];xx=np.column_stack([x,g[:,0]-g[:,0].mean()])
        enlarged=np.r_[best['theta'][:p],0.,best['theta'][p:2*p],0.,best['theta'][2*p:]]
        hh,_,_=derivatives(enlarged,xx,ww,y,d)
        bidx=[p,2*p+1];eidx=[i for i in range(len(enlarged)) if i not in bidx]
        he=hh[np.ix_(bidx,bidx)]-hh[np.ix_(bidx,eidx)]@np.linalg.solve(hh[np.ix_(eidx,eidx)],hh[np.ix_(eidx,bidx)])
        uc=components(enlarged,xx,ww,y,d)['scores'].sum(0)[bidx]
        stat=float(uc@np.linalg.solve(he,uc))
        score_differences.append(abs(stat-p0[0]['model_stat']))
        assert score_differences[-1]<2e-5
        # Independence model agrees with the sum of two correctly scaled score tests.
        fits.append(dict(columns=ww.shape[1],gradient=best['gradient'],nll=best['nll'],seconds=best['seconds']))
        fixture_checks.append(dict(columns=ww.shape[1],theta=best['theta'].tolist(),
             nll=best['nll'],model_p=p0[0]['model_p'],robust_p=p0[0]['robust_p']))
    assert max(flip_errors)<1e-10
    # Effects and CI code is exercised regardless of statistical significance.
    ff,est,se,model_se=full_effects(x,w,y,d,g[:,0])
    assert np.isfinite(est).all() and np.all(se>0) and np.all(model_se>0)
    result=dict(passed=True,engine_sha256=hashlib.sha256((ROOT/'src/model.py').read_bytes()).hexdigest(),gradient_max_error=gradient_error,
                independence_loglik_max_error=independent_error,
                hessian_max_error=max(h_errors),score_hessian_max_error=max(score_differences),
                allele_flip_max_error=max(flip_errors),null_fits=fits,
                elapsed_seconds=time.perf_counter()-began)
    (ROOT/'results/numerical_controls.json').write_text(json.dumps(result,indent=2)+'\n')
    # Public synthetic fixture for an independent base-R implementation.
    np.savetxt(ROOT/'results/synthetic_fixture.csv',np.column_stack([x[:,1:],y,d,g]),delimiter=',',
               header='age,exposure,continuous,binary,'+','.join('g'+str(j) for j in range(g.shape[1])),comments='')
    (ROOT/'results/fixture_python.json').write_text(json.dumps(dict(null_checks=fixture_checks,theta=ff['theta'].tolist(),
        nll=ff['nll'],effects=est.tolist(),sandwich_se=se.tolist(),model_se=model_se.tolist()),indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':run()
