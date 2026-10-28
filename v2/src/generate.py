"""Version two synthetic generator; nuisance intercept fixed across alternatives."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
import numpy as np
from simulate import intercept

def generate(seed,cell,rep,beta,stream=0):
    rng=np.random.default_rng(np.random.SeedSequence([seed,cell['id'],stream,rep]))
    n=cell['n'];maf=cell['maf']
    age=np.clip(rng.normal(size=n),-2.5,2.5);exposure=rng.binomial(1,.5,n)-.5
    x=np.column_stack([np.ones(n),age,exposure]);w=x[:,:2]
    g=rng.binomial(2,maf,size=(n,1)).astype(float)
    alpha=np.array([np.arctanh(.4),0.]) if cell['environment']=='constant' else np.array([np.arctanh(.3),.6])
    rho=np.tanh(w@alpha);z1=rng.normal(size=n);z2=rho*z1+np.sqrt(1-rho*rho)*rng.normal(size=n)
    a=intercept(cell['prevalence'],maf,0.)
    y=x@np.array([0.,.35,.25])+beta[0]*(g[:,0]-2*maf)+z1
    d=(x@np.array([a,.25,.4])+beta[1]*(g[:,0]-2*maf)+z2>0).astype(float)
    return dict(x=x,w=w,y=y,d=d,g=g,alpha=alpha,rho=rho,effects=beta,binary_intercept=a,maf=maf)
