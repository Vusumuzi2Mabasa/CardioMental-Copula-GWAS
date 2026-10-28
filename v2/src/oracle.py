"""Restricted nuisance score calculations; oracle information is never fitted.

The original engine remains unchanged. Active-index differentiation removes
known dependence coefficients from the nuisance projection, rather than merely
plugging known coefficients into a test that still estimates their uncertainty.
"""
from pathlib import Path
import sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2
from model import components,positive_inverse,fit,score_scan
from simulate import intercept

PRACTICAL=['Marginal covariance-aware','Constant model','Constant sandwich','Varying model','Varying sandwich']
ORACLES=['Known dependence model','Known dependence sandwich','All known model','All known sandwich']
METHODS=PRACTICAL+ORACLES

def true_theta(data):
    return np.r_[0.,.35,.25,data['binary_intercept'],.25,.4,0.,data['alpha']]

def restricted_fit(data):
    x,w,y,d=[data[k] for k in ['x','w','y','d']]
    if len(np.unique(d))!=2:raise ValueError('single binary class')
    # Marginal null fit supplies a deterministic start; alpha is held known.
    beta=np.linalg.lstsq(x,y,rcond=None)[0]
    from scipy.special import ndtri
    start=np.r_[beta,ndtri(d.mean()),0.,0.,np.log(np.std(y-x@beta))]
    def obj(v):
        c=components(np.r_[v,data['alpha']],x,w,y,d)
        return -c['loglik'].mean(),-c['scores'][:,:7].mean(0)
    bounds=[(-6,6)]*6+[(-2,2)]
    rr=minimize(obj,start,jac=True,method='L-BFGS-B',bounds=bounds,
                options=dict(ftol=1e-12,gtol=1e-7,maxiter=350,maxls=40))
    if not np.isfinite(rr.fun) or np.max(abs(rr.jac))>2e-5 or any(min(abs(v-a),abs(v-b))<1e-5 for v,(a,b) in zip(rr.x,bounds)):
        raise ValueError('restricted optimizer failed')
    return np.r_[rr.x,data['alpha']]

def restricted_score(theta,data,known_all=False):
    x,w,y,d=[data[k] for k in ['x','w','y','d']]
    c=components(theta,x,w,y,d)
    # Fully known oracle uses the known population genotype mean. Other tests
    # center at the sample mean as in the original engine with fitted intercepts.
    g=data['g'][:,0]-(2*data['maf'] if known_all else data['g'][:,0].mean())
    uindividual=g[:,None]*c['base']
    information=np.einsum('n,nkl->kl',g*g,c['h_base'])
    u=uindividual.sum(0)
    if not known_all:
        active=7; ds=np.empty((len(y),active,active));db=np.empty((len(y),2,active))
        for j in range(active):
            h=2e-5*(1+abs(theta[j]));a=theta.copy();b=theta.copy();a[j]+=h;b[j]-=h
            ca=components(a,x,w,y,d);cb=components(b,x,w,y,d)
            ds[:,:,j]=(ca['scores'][:,:active]-cb['scores'][:,:active])/(2*h)
            db[:,:,j]=(ca['base']-cb['base'])/(2*h)
        hh=-ds.sum(0);hh=(hh+hh.T)/2
        hbe=-np.einsum('n,nkp->kp',g,db);a=hbe@positive_inverse(hh)
        u-=a@c['scores'][:,:active].sum(0)
        uindividual-=c['scores'][:,:active]@a.T
        information-=a@hbe.T
    variance=uindividual.T@uindividual
    tm=float(u@positive_inverse(information)@u)
    tr=float(u@positive_inverse(variance)@u)
    return dict(model_p=float(chi2.sf(tm,2)),robust_p=float(chi2.sf(tr,2)),model_stat=tm,robust_stat=tr)

def evaluate(data):
    x,w,y,d,g=[data[k] for k in ['x','w','y','d','g']]
    out={}
    for label,ww in [('Marginal',w[:,:0]),('Constant',w[:,:1]),('Varying',w)]:
        start=time.perf_counter()
        targets=[('Marginal covariance-aware','robust_p')] if label=='Marginal' else [(label+' model','model_p'),(label+' sandwich','robust_p')]
        try:
            fitted=fit(x,ww,y,d);z=score_scan(fitted,x,ww,y,d,g)[0]
            if z['status']!='ok':raise ValueError(z['status'])
            for name,key in targets:out[name]=dict(p=z[key],status='ok',seconds=time.perf_counter()-start)
        except (ValueError,np.linalg.LinAlgError,FloatingPointError) as e:
            for name,key in targets:out[name]=dict(p=float('nan'),status=str(e),seconds=time.perf_counter()-start)
    for label,known_all in [('Known dependence',False),('All known',True)]:
        start=time.perf_counter()
        try:
            th=true_theta(data) if known_all else restricted_fit(data)
            z=restricted_score(th,data,known_all)
            for kind,key in [('model','model_p'),('sandwich','robust_p')]:
                out[label+' '+kind]=dict(p=z[key],status='ok',seconds=time.perf_counter()-start)
        except (ValueError,np.linalg.LinAlgError,FloatingPointError) as e:
            for kind in ['model','sandwich']:out[label+' '+kind]=dict(p=float('nan'),status=str(e),seconds=time.perf_counter()-start)
    return out
