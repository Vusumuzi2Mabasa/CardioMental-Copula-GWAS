"""Mixed phenotype marginal GWAS, explicit conjunction and single-causal coloc ABFs.

All inference uses observed phenotypes only. Generating truths are used solely
for evaluation. Whole-sample tagging evidence is not a shared-causal test.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import numpy as np
from scipy.special import expit, logsumexp, ndtr
from scipy.stats import norm, chi2

SCENARIOS=['null','cardio_only','mental_only','shared_same','shared_opposing','distinct_low','distinct_high','shared_plus_distinct']
NO_SHARED=['null','cardio_only','mental_only','distinct_low','distinct_high']

def generate(rng,hap,meta,scenario,regime,n=3200):
    # Independent whole haplotypes retain empirical local LD. No resampling of phenotypes.
    n_h,m=hap.shape;G=hap[rng.integers(n_h,size=n)]+hap[rng.integers(n_h,size=n)]
    age=np.clip(rng.normal(size=n),-2.5,2.5);sex=rng.binomial(1,.5,n)-.5
    exposure=rng.normal(size=n);X=np.column_stack([np.ones(n),age,sex,exposure,age*age])
    bc=np.zeros(m);bm=np.zeros(m);a,hi,lo=[meta[k] for k in ['anchor','high','low']]
    if scenario in ['cardio_only','shared_same','shared_opposing','distinct_low','distinct_high','shared_plus_distinct']:bc[a]=.18
    if scenario in ['mental_only','shared_same','shared_plus_distinct']:bm[a]=.35
    if scenario=='shared_opposing':bm[a]=-.35
    if scenario in ['distinct_low','shared_plus_distinct']:bm[lo]=.35
    if scenario=='distinct_high':bm[hi]=.35
    centered=G-2*hap.mean(0)
    z1=rng.normal(size=n);z2=rng.normal(size=n)
    rho=.35 if regime=='ideal' else np.tanh(np.arctanh(.30)+.35*age)
    e2=rho*z1+np.sqrt(1-rho*rho)*z2
    residual=z1 if regime=='ideal' else .70*z1+.45*rng.standard_t(5,size=n)/np.sqrt(5/3)+.15*(z1*z1-1)
    yc=.25*age+.15*sex+.20*exposure+.08*age*age+centered@bc+residual
    eta=-1.15+.20*age+.15*sex+.20*exposure+.05*age*age+centered@bm
    true_d=(eta+e2>0).astype(int);D=true_d.copy()
    treated=np.zeros(n,bool);keep=np.ones(n,bool)
    if regime=='realistic':
        # Treatment depends on blood pressure and age, not mental diagnosis.
        treated=rng.uniform(size=n)<expit(-1.8+.7*yc+.3*age)
        yc=yc-.6*treated+.25*rng.normal(size=n)
        sens,spec=.90,.98
        D=np.where(true_d,rng.uniform(size=n)<sens,rng.uniform(size=n)>spec).astype(int)
        # Joint measurement availability is MAR conditional on the fitted covariates.
        keep=rng.uniform(size=n)>expit(-2.7+.25*age+.20*exposure)
    # Clinical units are illustrative only; not estimated from an actual cohort.
    pressure=120+12*yc;cardio=(pressure>=140).astype(int)
    co=cardio*D
    odds=lambda c,d:(np.sum((c==1)&(d==1))+.5)*(np.sum((c==0)&(d==0))+.5)/((np.sum((c==1)&(d==0))+.5)*(np.sum((c==0)&(d==1))+.5))
    info=dict(n_generated=n,n_complete=int(keep.sum()),prevalence=float(D.mean()),true_prevalence=float(true_d.mean()),pressure_mean=float(pressure.mean()),pressure_sd=float(pressure.std()),treated_fraction=float(treated.mean()),missing_fraction=float(1-keep.mean()),comorbidity=float(co.mean()),comorbidity_or=float(odds(cardio,D)),true_shared=int(np.any((bc!=0)&(bm!=0))))
    return X[keep],yc[keep],D[keep],G[keep].astype(float),info,bc,bm

def linear(X,y,G):
    n,m=G.shape;p=X.shape[1]
    q=np.linalg.inv(X.T@X);gr=G-X@(q@(X.T@G));yr=y-X@(q@(X.T@y))
    den=np.sum(gr*gr,axis=0)
    if np.min(den)<1e-8:raise ValueError('monomorphic marker after adjustment')
    b=(gr.T@yr)/den;res=yr[:,None]-gr*b
    inf=gr*res/den
    robust=np.sum(inf*inf,axis=0)
    model=np.sum(res*res,axis=0)/(n-p-1)/den
    return dict(beta=b,var=robust,model_var=model,influence=inf)

def logistic(X,d,G):
    """Batched exact per-SNP logistic MLE and HC0 covariance, with explicit failures."""
    n,m=G.shape;p=X.shape[1];b=np.zeros((m,p+1));b[:,0]=np.log(d.mean()/(1-d.mean()))
    if d.sum()<5 or (1-d).sum()<5:raise ValueError('fewer than five cases/controls')
    for it in range(60):
        eta=X@b[:,:p].T+G*b[:,-1];mu=expit(eta);w=mu*(1-mu);e=d[:,None]-mu
        H=np.empty((m,p+1,p+1));H[:,:p,:p]=np.einsum('ni,nm,nj->mij',X,w,X,optimize=True)
        H[:,:p,p]=np.einsum('ni,nm->mi',X,w*G);H[:,p,:p]=H[:,:p,p]
        H[:,p,p]=np.sum(w*G*G,axis=0)
        score=np.column_stack([e.T@X,np.sum(e*G,axis=0)])
        delta=np.linalg.solve(H,score[...,None])[...,0]
        # Bounded Newton increments are only a solver safeguard; convergence still required.
        scale=np.maximum(1,np.max(np.abs(delta),axis=1)/2)
        b+=delta/scale[:,None]
        if np.max(np.abs(delta))<1e-8:break
    else:raise ValueError('logistic MLE did not converge')
    eta=X@b[:,:p].T+G*b[:,-1];mu=expit(eta);w=mu*(1-mu);e=d[:,None]-mu
    H[:,:p,:p]=np.einsum('ni,nm,nj->mij',X,w,X,optimize=True)
    H[:,:p,p]=np.einsum('ni,nm->mi',X,w*G);H[:,p,:p]=H[:,:p,p];H[:,p,p]=np.sum(w*G*G,axis=0)
    hi=np.linalg.inv(H);last=hi[:,-1,:]
    influence=e*(X@last[:,:p].T+G*last[:,-1])
    if np.max(np.abs(score))/n>1e-6:raise ValueError('logistic score tolerance')
    return dict(beta=b[:,-1],var=np.sum(influence**2,axis=0),model_var=hi[:,-1,-1],influence=influence,iterations=it+1)

def marginal(X,y,d,G):
    c=linear(X,y,G);b=logistic(X,d,G)
    z1=c['beta']/np.sqrt(c['var']);z2=b['beta']/np.sqrt(b['var'])
    p1=2*norm.sf(abs(z1));p2=2*norm.sf(abs(z2));cross=np.sum(c['influence']*b['influence'],axis=0)
    r=cross/np.sqrt(c['var']*b['var'])
    if np.max(abs(r))>=1:raise ValueError('singular cross-trait sandwich')
    joint=(z1*z1-2*r*z1*z2+z2*z2)/(1-r*r)
    return c,b,dict(p1=p1,p2=p2,iut=np.maximum(p1,p2),joint=chi2.sf(joint,2))

def labf(beta,var,sd):
    shrink=sd*sd/(sd*sd+var)
    return .5*(np.log1p(-shrink)+shrink*beta*beta/var)

def coloc(c,b,p12=1e-5):
    """Wakefield ABFs and standard coloc H0..H4 configuration enumeration.

    Continuous effect prior SD .15 (outcome scale), logistic .20 (log OR).
    Inputs are separate participant halves for primary inference.
    """
    a=labf(c['beta'],c['model_var'],.15);z=labf(b['beta'],b['model_var'],.20)
    s1,s2,both=logsumexp(a),logsumexp(z),logsumexp(a+z)
    # Exact off-diagonal enumeration avoids subtraction cancellation under near-identical signals.
    off=a[:,None]+z[None,:];np.fill_diagonal(off,-np.inf)
    logs=np.array([0,np.log(1e-4)+s1,np.log(1e-4)+s2,2*np.log(1e-4)+logsumexp(off),np.log(p12)+both])
    return np.exp(logs-logsumexp(logs))

def analyse_dataset(X,y,d,G,meta,copula=False):
    c,b,p=marginal(X,y,d,G);m=G.shape[1];a=meta['anchor'];n=len(y)
    # Fixed index partition independent of outcomes/genotypes; half sizes may differ by one.
    c1=linear(X[:n//2],y[:n//2],G[:n//2]);b2=logistic(X[n//2:],d[n//2:],G[n//2:])
    pp=coloc(c1,b2);overlap=coloc(c,b)
    row=dict(iut_p=float(min(1,m*p['iut'].min())),joint_p=float(min(1,m*p['joint'].min())),anchor_iut=float(p['iut'][a]),anchor_joint=float(p['joint'][a]),split_pp4=float(pp[4]),split_pp3=float(pp[3]),overlap_pp4=float(overlap[4]),prior_low_pp4=float(coloc(c1,b2,1e-6)[4]),prior_high_pp4=float(coloc(c1,b2,1e-4)[4]),min_p_cardio=float(p['p1'].min()),min_p_mental=float(p['p2'].min()))
    if copula:
        # Original v1/v2 Gaussian-probit engine, unchanged; diagnosis on realistic outcomes.
        import sys
        from pathlib import Path
        sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
        import importlib.util
        spec=importlib.util.spec_from_file_location('original_copula',Path(__file__).resolve().parents[2]/'src/model.py')
        orig=importlib.util.module_from_spec(spec);spec.loader.exec_module(orig)
        w=X[:,:2]
        try:
            fitted=orig.fit(X,w,y,d);scan=orig.score_scan(fitted,X,w,y,d,G[:,[a]])[0]
            if scan['status']!='ok':raise ValueError(scan['status'])
            row.update(copula_status='ok',copula_p=scan['robust_p'])
        except (ValueError,np.linalg.LinAlgError) as e:row.update(copula_status=str(e),copula_p=1.)
    else:row.update(copula_status='not_planned',copula_p=1.)
    fixture=dict(X=X,y=y,d=d,G=G,c_beta=c['beta'],c_var=c['model_var'],c_robust=c['var'],b_beta=b['beta'],b_var=b['model_var'],b_robust=b['var'],split_c_beta=c1['beta'],split_c_var=c1['model_var'],split_b_beta=b2['beta'],split_b_var=b2['model_var'],pp=pp)
    return row,fixture
