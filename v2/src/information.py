"""Post-protocol analytical diagnostic; never tunes the frozen main experiment.

Gaussian quadrature integrates the null Fisher information. The marginal
comparator uses Godambe information J.T V^-1 J for its stacked scores. This
calculation predicts asymptotic LOCAL power, not realised calibrated power.
"""
import json
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.special import ndtr
from scipy.stats import chi2,ncx2
from oracle import components
from generate import intercept
from run import ROOT,config,write

def information(environment,order=80,age_order=80):
    nodes,weights=hermgauss(order);weights=weights/np.sqrt(np.pi);nodes=np.sqrt(2)*nodes
    anodes,aweights=hermgauss(age_order);aweights=aweights/np.sqrt(np.pi)
    age=np.clip(np.sqrt(2)*anodes,-2.5,2.5)
    aa,ee,zz,dd=np.meshgrid(age,[-.5,.5],nodes,[0.,1.],indexing='ij')
    a=aa.ravel();e=ee.ravel();z=zz.ravel();d=dd.ravel()
    x=np.column_stack([np.ones(len(a)),a,e]);w=x[:,:2]
    alpha=np.array([np.arctanh(.4),0.]) if environment=='constant' else np.array([np.arctanh(.3),.6]) if environment=='varying' else np.zeros(2)
    b0=intercept(.15,.25,0.);theta=np.r_[0.,.35,.25,b0,.25,.4,0.,alpha]
    eta=x@theta[3:6];rho=np.tanh(w@alpha)
    probability=ndtr((2*d-1)*(eta+rho*z)/np.sqrt(1-rho*rho))
    mass=np.broadcast_to(aweights[:,None,None,None]*.5*weights[None,None,:,None],aa.shape).ravel()*probability
    c=components(theta,x,w,z+x@theta[:3],d)
    fisher=np.einsum('n,nij->ij',mass,c['h_base'])
    opg=(c['base']*mass[:,None]).T@c['base']
    scoremean=(mass[:,None]*c['base']).sum(0)
    # Marginal probit Fisher weight, averaged over covariates (not outcomes).
    ax,ex=np.meshgrid(age,[-.5,.5],indexing='ij');eta0=b0+.25*ax+.4*ex
    p=ndtr(eta0);phi=np.exp(-eta0**2/2)/np.sqrt(2*np.pi)
    b=phi**2/(p*(1-p));r=np.full_like(ax,.4) if environment=='constant' else np.tanh(np.arctanh(.3)+.6*ax) if environment=='varying' else np.zeros_like(ax)
    eb=np.sum(aweights[:,None]*.5*b);erb=np.sum(aweights[:,None]*.5*r*b)
    variance=np.array([[1.,erb],[erb,eb]]);sensitivity=np.diag([1.,eb])
    godambe=sensitivity@np.linalg.solve(variance,sensitivity)
    q=2*.25*.75
    return dict(joint=fisher*q,marginal=godambe*q,information_identity_error=float(np.max(abs(fisher-opg))),
                mean_score_max=float(max(abs(scoremean))),mass=float(mass.sum()),binary_intercept=b0)

def main():
    cfg=config();rows=[];controls=[];payload={}
    for env in ['independent','constant','varying']:
        a=information(env,80);b=information(env,120)
        agecheck=information(env,120,120)
        # Age clipping creates non-smooth quadrature, so report its sensitivity
        # separately from convergence over the smooth residual integral.
        difference=float(np.max(abs(a['joint']-b['joint'])))
        agedifference=float(np.max(abs(b['joint']-agecheck['joint'])))
        assert difference<2e-6,difference
        assert a['information_identity_error']<2e-5,a
        assert a['mean_score_max']<2e-5 and abs(a['mass']-1)<1e-12
        eig=np.linalg.eigvalsh(a['joint']-a['marginal'])
        assert eig[0]>-2e-7,eig
        if env=='independent':assert np.max(abs(a['joint']-a['marginal']))<1e-10
        controls.append(dict(environment=env,information_identity_error=a['information_identity_error'],mean_score_max=a['mean_score_max'],
            information_dominance_min_eigenvalue=float(eig[0]),residual_order80_vs120_max_difference=difference,
            age_order80_vs120_max_difference=agedifference))
        payload[env]=dict(joint=a['joint'].tolist(),marginal=a['marginal'].tolist(),binary_intercept=a['binary_intercept'])
        if env=='independent':continue
        for level,scale in zip(['weak','moderate'],cfg['scales']):
            for arch,base in cfg['architectures'].items():
                delta=np.array(base)*scale*np.sqrt(800)
                lj=float(delta@a['joint']@delta);lm=float(delta@a['marginal']@delta)
                rows.append(dict(environment=env,scenario=level+'_'+arch,joint_noncentrality=lj,marginal_noncentrality=lm,
                    information_ratio=lj/lm,joint_asymptotic_power=float(ncx2.sf(chi2.isf(.05,2),2,lj)),
                    marginal_asymptotic_power=float(ncx2.sf(chi2.isf(.05,2),2,lm))))
    write(ROOT/'results/information_predictions.csv',rows)
    result=dict(passed=True,scope='Post-protocol analytical diagnostic, not used for selection or calibration; local asymptotic power under correct Gaussian specification',
        quadrature_order=80,controls=controls,information=payload)
    (ROOT/'results/information_controls.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
