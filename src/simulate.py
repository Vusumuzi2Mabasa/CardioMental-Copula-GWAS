"""Prespecified synthetic outcomes and common-SNP panels; no patient data."""
from functools import lru_cache
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.special import ndtr,ndtri
from scipy.optimize import brentq
from scipy.stats import t as student_t

EFFECTS = dict(null=(0.,0.),continuous_only=(.20,0.),
               binary_only=(0.,.35),shared=(.20,.35),opposing=(.20,-.35))
ENVIRONMENTS = ['independent','constant','varying','t_varying']


@lru_cache(None)
def intercept(prevalence,maf,binary_effect):
    # Population quadrature, independent of every evaluation dataset.
    nodes,weights = hermgauss(60)
    age = np.clip(np.sqrt(2)*nodes,-2.5,2.5)
    weights = weights/np.sqrt(np.pi)
    shift=[]; mass=[]
    for exposure in [-.5,.5]:
        for g,pg in enumerate([(1-maf)**2,2*maf*(1-maf),maf*maf]):
            shift.extend(.25*age+.4*exposure+binary_effect*(g-2*maf))
            mass.extend(weights*.5*pg)
    shift=np.array(shift);mass=np.array(mass)
    return brentq(lambda a:np.sum(mass*ndtr(a+shift))-prevalence,-7,3)


def generate(seed,cell_id,rep,scenario,cell,markers=8):
    sid=list(EFFECTS).index(scenario)
    rng=np.random.default_rng(np.random.SeedSequence([seed,cell_id,rep,sid]))
    n=cell['n'];maf=cell['maf']; env=cell['environment']
    age=np.clip(rng.normal(size=n),-2.5,2.5)
    exposure=rng.binomial(1,.5,n)-.5
    x=np.column_stack([np.ones(n),age,exposure])
    w=np.column_stack([np.ones(n),age])
    if cell.get('ld',0):
        r=cell['ld']**np.abs(np.arange(markers)[:,None]-np.arange(markers))
        z=rng.normal(size=(n,2,markers))@np.linalg.cholesky(r).T
        g=(z<ndtri(maf)).sum(axis=1).astype(float)
    else:
        g=rng.binomial(2,maf,size=(n,markers)).astype(float)
    if env=='independent':alpha=np.array([0.,0.])
    elif env=='constant':alpha=np.array([np.arctanh(.4),0.])
    else:alpha=np.array([np.arctanh(.3),.6])
    rho=np.tanh(w@alpha)
    z1=rng.normal(size=n); z2=rho*z1+np.sqrt(1-rho*rho)*rng.normal(size=n)
    if env=='t_varying':
        scale=np.sqrt(rng.chisquare(4,n)/4)
        # Preserve normal margins while changing only the copula family.
        z1=ndtri(np.clip(student_t.cdf(z1/scale,4),1e-12,1-1e-12))
        z2=ndtri(np.clip(student_t.cdf(z2/scale,4),1e-12,1-1e-12))
    b1,b2=EFFECTS[scenario]
    continuous=.35*age+.25*exposure+b1*(g[:,0]-2*maf)+z1
    a=intercept(cell['prevalence'],maf,b2)
    latent=a+.25*age+.4*exposure+b2*(g[:,0]-2*maf)+z2
    binary=(latent>0).astype(float)
    return dict(x=x,w=w,y=continuous,d=binary,g=g,alpha=alpha,rho=rho,
                effects=(b1,b2),binary_intercept=a)
