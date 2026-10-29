"""Bivariate REML in the eigenspace of a projected relationship kernel.

Exact for two Gaussian traits under the specified mixed model. Gaussian working
likelihood for continuous/binary data; this approximation is evaluated separately.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

def matrices(t):
    a,b,c,d,e,f=t
    g11,g22=np.exp(2*a),np.exp(2*b);g12=np.exp(a+b)*np.tanh(c)
    r11,r22=np.exp(2*d),np.exp(2*e);r12=np.exp(d+e)*np.tanh(f)
    return np.array([[g11,g12],[g12,g22]]),np.array([[r11,r12],[r12,r22]])
def objective(t,lam,y):
    G,R=matrices(t);v1=lam*G[0,0]+R[0,0];v2=lam*G[1,1]+R[1,1];v12=lam*G[0,1]+R[0,1]
    det=v1*v2-v12*v12
    if np.any(det<=0):return np.inf,np.full(6,np.nan)
    i11=v2/det;i22=v1/det;i12=-v12/det
    q1=i11*y[:,0]+i12*y[:,1];q2=i12*y[:,0]+i22*y[:,1]
    f=.5*np.mean(np.log(det)+y[:,0]*q1+y[:,1]*q2)
    a11=i11-q1*q1;a22=i22-q2*q2;a12=i12-q1*q2
    dg=np.exp(t[0]+t[1])*(1-np.tanh(t[2])**2);dr=np.exp(t[3]+t[4])*(1-np.tanh(t[5])**2)
    grad=np.array([np.mean(lam*(a11*G[0,0]+a12*G[0,1])),np.mean(lam*(a22*G[1,1]+a12*G[0,1])),np.mean(lam*a12*dg),np.mean(a11*R[0,0]+a12*R[0,1]),np.mean(a22*R[1,1]+a12*R[0,1]),np.mean(a12*dr)])
    return f,grad
def starting(lam,y):
    A=np.array([[lam@lam,lam.sum()],[lam.sum(),len(lam)]])
    v=[]
    for j,k in [(0,0),(1,1),(0,1)]:v.append(np.linalg.solve(A,[np.sum(lam*y[:,j]*y[:,k]),np.sum(y[:,j]*y[:,k])]))
    total=np.mean(y*y,axis=0);g=np.maximum([v[0][0],v[1][0]],.1*total);r=np.maximum([v[0][1],v[1][1]],.1*total)
    gc=np.clip(v[2][0]/np.sqrt(g.prod()),-.8,.8);rc=np.clip(v[2][1]/np.sqrt(r.prod()),-.8,.8)
    return np.r_[.5*np.log(g),np.arctanh(gc),.5*np.log(r),np.arctanh(rc)]
def optimise(lam,y,start,null=False):
    ix=[0,1,3,4,5] if null else list(range(6));base=start.copy()
    if null:base[2]=0
    bounds=[(-7,2) if i not in [2,5] else (-3,3) for i in ix]
    def fn(v):
        t=base.copy();t[ix]=v;f,g=objective(t,lam,y);return f,g[ix]
    fit=minimize(fn,base[ix],jac=True,method='L-BFGS-B',bounds=bounds,options=dict(maxiter=300,ftol=1e-12,gtol=1e-7,maxls=40))
    t=base.copy();t[ix]=fit.x
    projected=fit.jac.copy()
    for j,(lo,hi) in enumerate(bounds):
        if fit.x[j]<=lo+1e-6 and projected[j]>0:projected[j]=0
        if fit.x[j]>=hi-1e-6 and projected[j]<0:projected[j]=0
    if not np.isfinite(fit.fun) or max(abs(projected))>2e-5:raise ValueError('REML convergence tolerance')
    return t,float(fit.fun)*len(lam),max(abs(projected))
def fit(lam,y,interval=True):
    start=starting(lam,y);tn,fn,gn=optimise(lam,y,start,True)
    t,fa,ga=optimise(lam,y,start,False)
    if fa>fn+1e-7:
        t,fa,ga=optimise(lam,y,tn,False)
    if fa>fn+1e-6:raise ValueError('Alternative REML worse than null')
    G,R=matrices(t);LR=max(0,2*(fn-fa));rg=float(np.tanh(t[2]));lo,hi=-1.,1.;ci_status='not_planned'
    if interval:
        h=1e-4;H=np.empty((6,6))
        for j in range(6):
            up=t.copy();dn=t.copy();up[j]+=h;dn[j]-=h
            H[:,j]=(objective(up,lam,y)[1]-objective(dn,lam,y)[1])/(2*h)*len(lam)
        H=(H+H.T)/2
        if np.linalg.eigvalsh(H).min()>1e-7:
            se=np.sqrt(np.linalg.inv(H)[2,2]);lo,hi=np.tanh(t[2]+np.array([-1,1])*1.95996398454*se);ci_status='ok'
        else:ci_status='nonpositive_information'
    result=dict(status='ok',genetic_c=float(G[0,0]),genetic_m=float(G[1,1]),genetic_covariance=float(G[0,1]),genetic_correlation=rg,residual_c=float(R[0,0]),residual_m=float(R[1,1]),residual_covariance=float(R[0,1]),residual_correlation=float(np.tanh(t[5])),LR=LR,asymptotic_p=float(chi2.sf(LR,1)),lower=float(lo),upper=float(hi),ci_status=ci_status,boundary=int(any(abs(t[j])>2.999 for j in [2,5]) or any(t[j]<-6.999 for j in [0,1,3,4])),max_gradient=float(max(ga,gn)),loglik_alt=-fa,loglik_null=-fn)
    return result,t,tn
