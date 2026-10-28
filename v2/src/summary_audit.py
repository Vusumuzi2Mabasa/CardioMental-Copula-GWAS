"""Independently reconstruct reported decisions and audit bootstrap summaries."""
import json
import numpy as np
from scipy.stats import binomtest,chi2,ncx2
from run import ROOT,config,cells,scenarios,read
from data import load
from oracle import METHODS

def close(a,b,tol=1e-12):
    assert np.isclose(float(a),float(b),atol=tol,rtol=0,equal_nan=True),(a,b)

def main():
    cfg=config();summary=read(ROOT/'results/summary.csv');pairs=read(ROOT/'results/paired_contrasts.csv')
    intervals=read(ROOT/'results/calibration_aware_intervals.csv')
    threshold=json.loads((ROOT/'results/thresholds.json').read_text())['thresholds']
    draws=np.load(ROOT/'results/bootstrap_draws.npz',allow_pickle=False)
    checked=0;checked_pairs=0;replayed=0
    for cell in cells():
        ss=scenarios(cell);cal=load(cell,'calibration','null');aa={s['name']:load(cell,'evaluation',s['name']) for s in ss}
        boot=draws['cell_'+str(cell['id'])];assert boot.shape==(cfg['bootstrap_replicates'],9,9)
        rng=np.random.default_rng(np.random.SeedSequence([cfg['seed_uncertainty'],cell['id']]))
        # Replay RNG draws for all iterations, independently calculate selected
        # first and last draws with per-method sorted calibration samples.
        for b in range(cfg['bootstrap_replicates']):
            ix=rng.integers(0,cal['n'],cal['n'])
            keep=b in [0,1,cfg['bootstrap_replicates']-1]
            if keep:
                rank=int(.05*(cal['n']+1))-1
                cuts=np.array([sorted(cal['p'][ix,j])[rank] if cal['good'][:,j].all() else np.nan for j in range(9)])
            for sid,s in enumerate(ss):
                a=aa[s['name']];ex=rng.integers(0,a['n'],a['n'])
                if keep:
                    expected=np.array([np.mean((a['p'][ex,j]<cuts[j])&a['good'][ex,j]) if np.isfinite(cuts[j]) else np.nan for j in range(9)])
                    assert np.allclose(expected,boot[b,sid],atol=1e-12,rtol=0,equal_nan=True);replayed+=1
        for sid,s in enumerate(ss):
            a=aa[s['name']];n=a['n'];decisions={}
            for alpha in cfg['alpha_grid']:
                for j,m in enumerate(METHODS):
                    r=next(z for z in summary if int(z['id'])==cell['id'] and z['scenario']==s['name'] and z['method']==m and float(z['alpha'])==alpha)
                    t=next(z for z in threshold if z['cell']==cell['id'] and z['method']==m and z['alpha']==alpha)
                    cut=t['cutoff'];nom=(a['p'][:,j]<alpha)&a['good'][:,j]
                    dec=(a['p'][:,j]<cut)&a['good'][:,j] if cut is not None else None
                    assert int(r['attempted'])==n and int(r['valid'])==a['good'][:,j].sum() and int(r['failures'])==int((~a['good'][:,j]).sum())
                    close(r['nominal_rate'],nom.mean())
                    if dec is not None:
                        close(r['calibrated_rate'],dec.mean())
                        ci=binomtest(int(dec.sum()),n).proportion_ci(method='wilson')
                        close(r['conditional_lower'],ci.low);close(r['conditional_upper'],ci.high)
                        close(r['conditional_MCSE'],np.sqrt(dec.mean()*(1-dec.mean())/n))
                    if alpha==.05:
                        decisions[m]=dec
                        bb=next(z for z in intervals if int(z['id'])==cell['id'] and z['scenario']==s['name'] and z['method']==m)
                        q=np.quantile(boot[:,sid,j],[.025,.975]);close(bb['lower'],q[0]);close(bb['upper'],q[1])
                    checked+=1
            for row in [z for z in pairs if int(z['id'])==cell['id'] and z['scenario']==s['name']]:
                j=METHODS.index(row['method']);k=METHODS.index(row['reference'])
                if decisions[row['method']] is not None and decisions[row['reference']] is not None:
                    diff=decisions[row['method']].astype(int)-decisions[row['reference']].astype(int)
                    close(row['difference'],diff.mean());close(row['conditional_MCSE'],diff.std(ddof=1)/np.sqrt(n))
                q=np.quantile(boot[:,sid,j]-boot[:,sid,k],[.025,.975])
                close(row['calibration_aware_lower'],q[0]);close(row['calibration_aware_upper'],q[1]);checked_pairs+=1
    assert checked==len(summary)==1620 and checked_pairs==len(pairs)==216
    matrices=json.loads((ROOT/'results/information_controls.json').read_text())['information']
    for row in read(ROOT/'results/information_predictions.csv'):
        sc=next(s for s in scenarios(cells()[0]) if s['name']==row['scenario']);delta=np.array(sc['beta'])*np.sqrt(800)
        for key,matrix in [('joint','joint'),('marginal','marginal')]:
            lam=float(delta@np.array(matrices[row['environment']][matrix])@delta)
            close(row[key+'_noncentrality'],lam)
            close(row[key+'_asymptotic_power'],ncx2.sf(chi2.isf(.05,2),2,lam))
    result=dict(passed=True,summary_rows=checked,paired_rows=checked_pairs,replayed_bootstrap_scenario_draws=replayed,
                Wilson_checked_with_independent_scipy=True,all_bootstrap_quantiles_verified=True,information_predictions_verified=True)
    (ROOT/'results/summary_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
if __name__=='__main__':main()
