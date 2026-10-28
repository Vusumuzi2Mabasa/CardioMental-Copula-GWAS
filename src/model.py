"""Exact Gaussian/probit joint likelihood and nuisance-adjusted SNP score tests.

No covariance matrix, copula probability or test statistic is clipped to force
success. Information failures are returned explicitly to the simulation runner.
"""
import os
for _k in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_k] = '1'
import time
import numpy as np
from scipy.special import log_ndtr, ndtr, ndtri, erfcx
from scipy.optimize import minimize
from scipy.stats import chi2

LOG2PI = np.log(2*np.pi)
METHODS = ['Marginal Bonferroni', 'Marginal covariance-aware',
           'Independence likelihood', 'Constant model score',
           'Constant sandwich score', 'Varying model score',
           'Varying sandwich score']


def components(theta, x, w, y, d):
    p = x.shape[1]
    mu = x @ theta[:p]
    eta = x @ theta[p:2*p]
    sigma = np.exp(theta[2*p])
    rho = np.tanh(w @ theta[2*p+1:]) if w.shape[1] else np.zeros(len(y))
    s = np.sqrt(1-rho*rho)
    z = (y-mu)/sigma
    a = (eta+rho*z)/s
    t = 2*d-1
    logprob = log_ndtr(t*a)
    r=t*a
    mills=np.empty_like(r)
    negative=r<0
    mills[negative]=np.sqrt(2/np.pi)/erfcx(-r[negative]/np.sqrt(2))
    mills[~negative]=np.exp(-.5*r[~negative]**2-.5*LOG2PI-logprob[~negative])
    k=t*mills
    loglik = -.5*LOG2PI-np.log(sigma)-.5*z*z+logprob
    u_mu = (z-k*rho/s)/sigma
    u_eta = k/s
    u_logsigma = z*z-1-k*rho*z/s
    u_alpha = (k*(z+rho*eta)/s)[:,None]*w
    scores = np.column_stack([u_mu[:,None]*x, u_eta[:,None]*x,
                              u_logsigma, u_alpha])
    v = mills*(r+mills)
    tail=r < -30
    u=-r[tail]
    delta=1/u-2/u**3+10/u**5-74/u**7
    v[tail]=(u+delta)*delta
    h_base = np.zeros((len(y),2,2))
    h_base[:,0,0] = (1+v*rho*rho/(s*s))/(sigma*sigma)
    h_base[:,0,1] = h_base[:,1,0] = -rho*v/(sigma*s*s)
    h_base[:,1,1] = v/(s*s)
    return dict(loglik=loglik, scores=scores,
                base=np.column_stack([u_mu,u_eta]), h_base=h_base,
                sigma=sigma, rho=rho, eta=eta, z=z)


def objective(theta, x, w, y, d):
    c = components(theta,x,w,y,d)
    return -np.mean(c['loglik']), -np.mean(c['scores'],axis=0)


def fit(x, w, y, d, multistart=False):
    started = time.perf_counter()
    if len(np.unique(d)) != 2:
        raise ValueError('binary outcome has fewer than two classes')
    p = x.shape[1]
    beta = np.linalg.lstsq(x,y,rcond=None)[0]
    sigma = np.sqrt(np.mean((y-x@beta)**2))
    gamma = np.zeros(p); gamma[0] = ndtri(np.mean(d))
    start = np.r_[beta,gamma,np.log(sigma),np.zeros(w.shape[1])]
    bounds = [(-6,6)]*(2*p)+[(-2,2)]+[(-2.5,2.5)]*w.shape[1]
    starts = [start]
    if multistart and w.shape[1]:
        for rho in [-.4,.4]:
            z = start.copy(); z[2*p+1] = np.arctanh(rho); starts.append(z)
    options = dict(ftol=1e-12,gtol=1e-7,maxiter=350,maxls=40)
    solutions = [minimize(objective,s,args=(x,w,y,d),jac=True,
                          method='L-BFGS-B',bounds=bounds,options=options)
                 for s in starts]
    best = min(solutions,key=lambda z:z.fun)
    # Deterministic rescue for a poorly converged first start.
    if np.max(np.abs(best.jac)) > 2e-5 and not multistart and w.shape[1]:
        rescue = start.copy(); rescue[2*p+1] = .4
        solutions.append(minimize(objective,rescue,args=(x,w,y,d),jac=True,
                         method='L-BFGS-B',bounds=bounds,options=options))
        best = min(solutions,key=lambda z:z.fun)
    gradient = float(np.max(np.abs(best.jac)))
    boundary = any(abs(v-b[0])<1e-5 or abs(v-b[1])<1e-5
                   for v,b in zip(best.x,bounds))
    if not np.isfinite(best.fun) or gradient > 2e-5 or boundary:
        raise ValueError(f'null/full optimizer: gradient={gradient:.4g}, boundary={boundary}')
    return dict(theta=best.x,nll=float(best.fun)*len(y),gradient=gradient,
                iterations=int(best.nit),starts=len(solutions),
                seconds=time.perf_counter()-started)


def derivatives(theta,x,w,y,d,step=2e-5):
    n = len(y); p = len(theta)
    ds = np.empty((n,p,p)); db = np.empty((n,2,p))
    for j in range(p):
        h = step*(1+abs(theta[j]))
        a = theta.copy(); b = theta.copy(); a[j]+=h; b[j]-=h
        ca = components(a,x,w,y,d); cb = components(b,x,w,y,d)
        ds[:,:,j] = (ca['scores']-cb['scores'])/(2*h)
        db[:,:,j] = (ca['base']-cb['base'])/(2*h)
    raw = -ds.sum(0)
    hessian = (raw+raw.T)/2
    return hessian, db, float(np.max(np.abs(raw-raw.T)))


def positive_inverse(a):
    a = (a+a.T)/2
    e = np.linalg.eigvalsh(a)
    if not np.isfinite(e).all() or e[0] <= max(1e-10,1e-10*e[-1]):
        raise ValueError('non-positive or ill-conditioned information')
    return np.linalg.inv(a)


def score_scan(fitted,x,w,y,d,genotypes):
    theta = fitted['theta']; c = components(theta,x,w,y,d)
    h,db,asym = derivatives(theta,x,w,y,d)
    hi = positive_inverse(h)
    nuisance_score = c['scores'].sum(0)
    output = []
    for raw_g in genotypes.T:
        g = raw_g-np.mean(raw_g)
        if np.var(g) == 0:
            output.append(dict(status='monomorphic')); continue
        try:
            hbe = -np.einsum('n,nkp->kp',g,db)
            adjustment = hbe@hi
            individual = g[:,None]*c['base']-c['scores']@adjustment.T
            u = np.sum(g[:,None]*c['base'],axis=0)-adjustment@nuisance_score
            information = np.einsum('n,nkl->kl',g*g,c['h_base'])-adjustment@hbe.T
            robust = individual.T@individual
            im = positive_inverse(information); ir = positive_inverse(robust)
            tm = float(u@im@u); tr = float(u@ir@u)
            marginal = u*u/np.diag(robust)
            output.append(dict(status='ok',model_p=float(chi2.sf(tm,2)),
                robust_p=float(chi2.sf(tr,2)),model_stat=tm,robust_stat=tr,
                bonf_p=float(min(1,2*np.min(chi2.sf(marginal,1)))),
                hessian_asymmetry=asym))
        except (ValueError,np.linalg.LinAlgError) as e:
            output.append(dict(status=str(e)))
    return output


def full_effects(x,w,y,d,g):
    xx = np.column_stack([x,g-np.mean(g)])
    fitted = fit(xx,w,y,d)
    c = components(fitted['theta'],xx,w,y,d)
    h,_,_ = derivatives(fitted['theta'],xx,w,y,d)
    hi = positive_inverse(h)
    sandwich = hi@(c['scores'].T@c['scores'])@hi
    p = xx.shape[1]
    ix = [p-1,2*p-1]
    return fitted, fitted['theta'][ix], np.sqrt(np.diag(sandwich)[ix]), np.sqrt(np.diag(hi)[ix])
