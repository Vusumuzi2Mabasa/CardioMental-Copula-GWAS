"""Independent gradient, dense likelihood and allele/kernel invariance controls."""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.optimize._numdiff import approx_derivative
from scipy.stats import multivariate_normal
from reml import objective,matrices,fit
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'v3/src'))
from prepare import dump
def main():
    rng=np.random.default_rng(29103199);lam=np.linspace(.05,3,80);t=np.array([-.4,-.6,.3,-.2,-.4,.2]);G,R=matrices(t)
    y=np.array([rng.multivariate_normal([0,0],v*G+R) for v in lam])
    f,g=objective(t,lam,y);gd=approx_derivative(lambda u:np.array([objective(u,lam,y)[0]]),t).ravel()
    error=max(abs(g-gd));assert error<1e-7
    direct=-sum(multivariate_normal.logpdf(yy,mean=[0,0],cov=v*G+R) for yy,v in zip(y,lam))/len(lam)-np.log(2*np.pi)
    assert abs(direct-f)<1e-12
    result,theta,tn=fit(lam,y);flipped,_,_=fit(lam,y*np.array([1,-1]))
    assert abs(result['genetic_correlation']+flipped['genetic_correlation'])<1e-8
    assert abs(result['LR']-flipped['LR'])<1e-8
    p=ROOT/'v4/results/tests';p.mkdir(parents=True,exist_ok=True)
    dump(p/'reml_fixture.json',dict(eigenvalues=lam.tolist(),y=y.tolist(),theta=theta.tolist(),theta_null=tn.tolist(),expected=result))
    out=dict(passed=True,gradient_max_error=float(error),dense_loglikelihood_error=float(abs(direct-f)),outcome_sign_flip=True)
    dump(p/'numerical_controls.json',out);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
