"""Evaluate genetic correlation estimates and independent empirical calibration."""
from pathlib import Path
import sys,csv,json,argparse
import numpy as np
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[2];sys.path.append(str(ROOT/'v3/src'))
from prepare import dump,table
from polygenic import tag,V
def load(phase,validation=False):
    rows=[]
    base=V/('results/polygenic_validation' if validation else 'results')
    for p in sorted((base/'polygenic_chunks').glob(phase+'_*.csv')):
        with p.open() as f:
            for r in csv.DictReader(f):
                for k,v in list(r.items()):
                    if k not in ['population','phase','mode','status','ci_status']:r[k]=float(v)
                rows.append(r)
    return rows
def cell(r):return tag(r['population'],int(r['n']),r['residual_truth'])+'_'+r['mode']
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--validation',action='store_true');args=ap.parse_args()
    base=V/('results/polygenic_validation' if args.validation else 'results')
    nc,ne=(1000,200) if args.validation else (200,100);rank=nc-int(np.floor(.05*(nc+1)))
    cal=load('cal',args.validation);rows=load('eval',args.validation);assert len(cal)==16*nc and len(rows)==48*ne
    ts=json.loads((base/'polygenic_thresholds.json').read_text())['thresholds'];summ=[]
    factor=.88*norm.pdf(norm.isf(.15))
    for key,t in ts.items():
        bank=[r for r in cal if cell(r)==key];assert len(bank)==nc
        assert t['LR_cutoff']==float(np.sort([r['LR'] for r in bank])[rank])
        for truth in [-.5,0.,.5]:
            a=[r for r in rows if cell(r)==key and r['genetic_truth']==truth];assert len(a)==ne
            valid=[r for r in a if r['status']=='ok'];ci=[r for r in valid if r['ci_status']=='ok'];r=a[0]
            scale=factor if r['mode']=='observed_binary' else 1.
            summ.append(dict(population=r['population'],n=int(r['n']),residual_truth=r['residual_truth'],mode=r['mode'],genetic_truth=truth,attempted=ne,valid=len(valid),failures=ne-len(valid),boundary=sum(r['boundary'] for r in valid),valid_intervals=len(ci),mean_rg=float(np.mean([r['genetic_correlation'] for r in valid])),bias_rg=float(np.mean([r['genetic_correlation']-truth for r in valid])),rmse_rg=float(np.sqrt(np.mean([(r['genetic_correlation']-truth)**2 for r in valid]))),mean_covariance_liability_approx=float(np.mean([r['genetic_covariance']/scale for r in valid])),mean_variance_m_liability_approx=float(np.mean([r['genetic_m']/scale**2 for r in valid])),coverage_attempted=sum(r['lower']<=truth<=r['upper'] for r in ci)/ne,coverage_valid=sum(r['lower']<=truth<=r['upper'] for r in ci)/len(ci) if ci else 0,nominal_rate=sum(r['status']=='ok' and r['asymptotic_p']<.05 for r in a)/ne,calibrated_rate=sum(t['active'] and r['status']=='ok' and r['LR']>t['LR_cutoff'] for r in a)/ne,calibration_active=t['active'],mean_phenotypic_correlation=float(np.mean([r['phenotypic_correlation' if r['mode']=='observed_binary' else 'liability_phenotypic_correlation'] for r in a]))))
    table(base/'polygenic_rates.csv',summ)
    # Paired bootstrap: both outcome modes share each generated cohort and replicate index.
    rng=np.random.default_rng(29103113);B=400;draws=[]
    cb={key:sorted([r for r in cal if cell(r)==key],key=lambda r:r['replicate']) for key in ts}
    eb={(key,g):sorted([r for r in rows if cell(r)==key and r['genetic_truth']==g],key=lambda r:r['replicate']) for key in ts for g in [-.5,0.,.5]}
    for bi in range(B):
        for n in [600,1200]:
            for truth in [-.5,0.,.5]:
                counts={'liability':[],'observed_binary':[]}
                for pop in ['YRI','LWK']:
                    for re in [0.,.4]:
                        ii=rng.integers(nc,size=nc);jj=rng.integers(ne,size=ne)
                        for mode in counts:
                            key=tag(pop,n,re)+'_'+mode
                            bank=cb[key];evals=eb[key,truth]
                            cut=np.sort([bank[i]['LR'] for i in ii])[rank]
                            counts[mode].append(np.mean([ts[key]['active'] and evals[j]['status']=='ok' and evals[j]['LR']>cut for j in jj]))
                for mode,values in counts.items():draws.append(dict(bootstrap=bi,n=n,genetic_truth=truth,mode=mode,rate=float(np.mean(values))))
    table(base/'polygenic_bootstrap.csv',draws)
    agg=[]
    for n in [600,1200]:
        for truth in [-.5,0.,.5]:
            for mode in ['liability','observed_binary']:
                a=[r for r in summ if r['n']==n and r['genetic_truth']==truth and r['mode']==mode]
                dd=[r['rate'] for r in draws if r['n']==n and r['genetic_truth']==truth and r['mode']==mode];lo,hi=np.quantile(dd,[.025,.975])
                agg.append(dict(n=n,genetic_truth=truth,mode=mode,attempted=ne*4,failures=sum(r['failures'] for r in a),**{k:float(np.mean([r[k] for r in a])) for k in ['mean_rg','bias_rg','rmse_rg','coverage_attempted','nominal_rate','calibrated_rate']},bootstrap_lower=float(lo),bootstrap_upper=float(hi)))
    table(base/'polygenic_aggregate.csv',agg)
    dump(base/'polygenic_audit.json',dict(passed=True,unique_calibration_cohorts=8*nc,unique_evaluation_cohorts=24*ne,outcome_fits=16*nc+48*ne,outcome_modes_not_independent=True,calibration_fit_failures=sum(r['status']!='ok' for r in cal),evaluation_fit_failures=sum(r['status']!='ok' for r in rows),inactive_calibration_strata=sum(not t['active'] for t in ts.values()),summary_rows=len(summ),aggregate_rows=len(agg),bootstrap_replicates=B,binary_likelihood='Gaussian working approximation; errors and coverage evaluated separately from exact liability benchmark'))
    print('Polygenic summaries complete',flush=True)
if __name__=='__main__':main()
