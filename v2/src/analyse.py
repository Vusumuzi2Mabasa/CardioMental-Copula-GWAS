"""Conditional and calibration-aware uncertainty from paired held-out outputs."""
import json,platform,sys
import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run import ROOT,config,cells,scenarios,write,read
from data import load
from oracle import METHODS

SHORT=['Marginal','C model','C robust','V model','V robust','Known rho model','Known rho robust','All known model','All known robust']
SHOW=[0,2,4,6,8]
COLORS=['#222222','#4679a6','#29958c','#a05b9c','#bb4949']
plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':10})

def ci(p,n):
    z=norm.ppf(.975);den=1+z*z/n;c=(p+z*z/(2*n))/den;h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return float(max(0,c-h)),float(min(1,c+h))
def meta(cell):return {k:cell[k] for k in ['id','n','environment']}
def save(fig,name):
    fig.savefig(ROOT/'figures'/name,dpi=200,bbox_inches='tight',facecolor='white');plt.close(fig)
def main():
    audit=json.loads((ROOT/'results/audit.json').read_text());assert audit['passed']
    cfg=config();frozen=json.loads((ROOT/'results/thresholds.json').read_text())['thresholds']
    summaries=[];pairs=[];bootrows=[];banks=[];runtime=[];lookup={};bootlookup={}
    contrasts=[(4,0),(4,2),(6,4),(8,4),(7,3),(5,3)]
    for cell in cells():
        cal=load(cell,'calibration','null');ss=scenarios(cell)
        arrays={s['name']:load(cell,'evaluation',s['name']) for s in ss}
        rng=np.random.default_rng(np.random.SeedSequence([cfg['seed_uncertainty'],cell['id']]))
        nboot=cfg['bootstrap_replicates'];br={s['name']:np.full((nboot,len(METHODS)),np.nan) for s in ss}
        # Resample calibration datasets jointly across methods, then evaluation
        # datasets jointly across methods. Separate streams remain independent.
        for b in range(nboot):
            ix=rng.integers(0,cal['n'],cal['n']);rank=int(np.floor(.05*(cal['n']+1)))
            cut=np.partition(cal['p'][ix],rank-1,axis=0)[rank-1]
            cut[~cal['good'].all(0)]=np.nan
            for s in ss:
                a=arrays[s['name']];ex=rng.integers(0,a['n'],a['n'])
                br[s['name']][b]=np.mean((a['p'][ex]<cut)&a['good'][ex],axis=0)
                br[s['name']][b,~np.isfinite(cut)]=np.nan
        for s in ss:
            name=s['name'];a=arrays[name];n=a['n'];valid=a['good'].sum(0);fail=n-valid
            for alpha in cfg['alpha_grid']:
                cut=np.array([next(t['cutoff'] for t in frozen if t['cell']==cell['id'] and t['method']==m and t['alpha']==alpha) for m in METHODS],dtype=float)
                nom=(a['p']<alpha)&a['good'];dec=(a['p']<cut)&a['good']
                for j,m in enumerate(METHODS):
                    rate=float(dec[:,j].mean()) if np.isfinite(cut[j]) else float('nan');nomrate=float(nom[:,j].mean())
                    lo,hi=ci(rate,n) if np.isfinite(rate) else (float('nan'),float('nan'))
                    row=dict(**meta(cell),scenario=name,method=m,alpha=alpha,attempted=n,valid=int(valid[j]),failures=int(fail[j]),
                        cutoff=cut[j],nominal_rate=nomrate,calibrated_rate=rate,
                        valid_nominal_rate=float(nom[:,j].sum()/valid[j]) if valid[j] else float('nan'),
                        valid_calibrated_rate=float(dec[:,j].sum()/valid[j]) if valid[j] and np.isfinite(cut[j]) else float('nan'),
                        conditional_MCSE=float(np.sqrt(rate*(1-rate)/n)),
                        conditional_lower=lo,conditional_upper=hi,nominal_lower=ci(nomrate,n)[0],nominal_upper=ci(nomrate,n)[1])
                    summaries.append(row)
                    if alpha==.05:lookup[cell['id'],name,j]=row
                if alpha==.05:
                    for j,k in contrasts:
                        d=dec[:,j].astype(float)-dec[:,k];estimate=float(d.mean());se=float(d.std(ddof=1)/np.sqrt(n))
                        if not np.isfinite(cut[[j,k]]).all():estimate=se=float('nan')
                        quant=np.quantile(br[name][:,j]-br[name][:,k],[.025,.975])
                        pairs.append(dict(**meta(cell),scenario=name,method=METHODS[j],reference=METHODS[k],difference=estimate,
                            conditional_MCSE=se,conditional_lower=estimate-1.96*se,conditional_upper=estimate+1.96*se,
                            calibration_aware_lower=quant[0],calibration_aware_upper=quant[1],bootstrap_replicates=nboot))
            for j,m in enumerate(METHODS):
                lo,hi=np.quantile(br[name][:,j],[.025,.975]);r=dict(**meta(cell),scenario=name,method=m,
                    lower=lo,upper=hi,bootstrap_SD=float(br[name][:,j].std(ddof=1)),replicates=nboot)
                bootrows.append(r);bootlookup[cell['id'],name,j]=r
            for bank in range(cfg['bank_partitions']):
                ix=np.arange(bank,cal['n'],cfg['bank_partitions']);rank=int(np.floor(.05*(len(ix)+1)))
                cuts=np.partition(cal['p'][ix],rank-1,axis=0)[rank-1]
                for j,m in enumerate(METHODS):
                    failures=int((~cal['good'][ix,j]).sum())
                    banks.append(dict(**meta(cell),scenario=name,bank=bank,calibration_n=len(ix),method=m,
                        cutoff=cuts[j] if not failures else float('nan'),calibration_failures=failures,
                        evaluation_rate=float(np.mean((a['p'][:,j]<cuts[j])&a['good'][:,j])) if not failures else float('nan')))
        for j,m in enumerate(METHODS):runtime.append(dict(**meta(cell),method=m,mean_seconds=float(arrays['null']['seconds'][:,j].mean())))
        print('Summarised cell',cell['id'],flush=True)
    for name,rows in [('summary.csv',summaries),('paired_contrasts.csv',pairs),('calibration_aware_intervals.csv',bootrows),('subbank_sensitivity.csv',banks),('runtime.csv',runtime)]:write(ROOT/'results'/name,rows)
    # Main nominal null calibration includes model and sandwich versions.
    fig,axes=plt.subplots(2,2,figsize=(13,8),sharey=True)
    for ax,cell in zip(axes.flat,cells()):
        for j in range(len(METHODS)):
            r=lookup[cell['id'],'null',j];ax.errorbar(j,r['nominal_rate'],yerr=[[r['nominal_rate']-r['nominal_lower']],[r['nominal_upper']-r['nominal_rate']]],fmt='o',color='#356e92',capsize=2)
        ax.axhline(.05,color='gray',ls='--');ax.set_xticks(range(9),SHORT,rotation=55,ha='right');ax.set_ylabel('Nominal null rejection')
        limits=[lookup[c['id'],'null',j] for c in cells() for j in range(len(METHODS))]
        ax.set_ylim(max(0,min(.025,min(r['nominal_lower'] for r in limits))-.003),max(.085,max(r['nominal_upper'] for r in limits))+.003)
        ax.set_title(f"n={cell['n']} | {cell['environment']}")
    fig.suptitle('SIMULATED data | 10,000 held-out nulls per setting | 95% Wilson intervals');fig.tight_layout();save(fig,'null_calibration.png')
    altnames=[s['name'] for s in scenarios(cells()[0])][1:]
    labels=['W C','W B','W S','W O','M C','M B','M S','M O']
    fig,axes=plt.subplots(2,2,figsize=(13,8),sharey=True)
    for ax,cell in zip(axes.flat,cells()):
        for j,color in zip(SHOW,COLORS):
            r=[lookup[cell['id'],s,j] for s in altnames]
            xx=np.arange(8)+(SHOW.index(j)-2)*.09
            yy=[z['calibrated_rate'] for z in r]
            ax.errorbar(xx,yy,yerr=[[z['calibrated_rate']-z['conditional_lower'] for z in r],[z['conditional_upper']-z['calibrated_rate'] for z in r]],fmt='o',ms=3,capsize=1,color=color,label=SHORT[j])
        ax.set_xticks(range(8),labels,rotation=25);ax.set_ylabel('Calibrated association power');ax.set_ylim(0,1.03);ax.set_title(f"n={cell['n']} | {cell['environment']}")
    axes[0,0].legend(ncol=2,loc='upper left');fig.suptitle('SIMULATED local alternatives | W weak M moderate | C continuous B binary S shared O opposing')
    fig.tight_layout();save(fig,'moderate_power.png')
    fig,axes=plt.subplots(2,2,figsize=(13,8),sharey=True)
    for ax,cell in zip(axes.flat,cells()):
        for (j,k),color in zip([(4,0),(4,2),(6,4),(8,4)],COLORS):
            rr=[next(r for r in pairs if r['id']==cell['id'] and r['scenario']==s and r['method']==METHODS[j] and r['reference']==METHODS[k]) for s in altnames]
            xx=np.arange(8)+([(4,0),(4,2),(6,4),(8,4)].index((j,k))-1.5)*.12
            # Percentile intervals can fail to contain the point estimate; draw
            # endpoints directly instead of forcing nonnegative error lengths.
            ax.vlines(xx,[r['calibration_aware_lower'] for r in rr],[r['calibration_aware_upper'] for r in rr],color=color,lw=1)
            ax.scatter(xx,[r['difference'] for r in rr],s=12,color=color,label=SHORT[j]+' minus '+SHORT[k])
        ax.axhline(0,color='gray',ls='--');ax.set_xticks(range(8),labels,rotation=25);ax.set_ylabel('Paired power difference');ax.set_title(f"n={cell['n']} | {cell['environment']}")
    axes[0,0].legend(fontsize=9);fig.suptitle('SIMULATED data | 95% paired bootstrap intervals include calibration and evaluation sampling')
    fig.tight_layout();save(fig,'oracle_gaps.png')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),sharey=True)
    for ax,cell in zip(axes,[cells()[1],cells()[3]]):
        for j,color in zip(SHOW,COLORS):
            xx=[];yy=[]
            for alpha in cfg['alpha_grid']:
                null=next(r for r in summaries if r['id']==cell['id'] and r['method']==METHODS[j] and r['scenario']=='null' and r['alpha']==alpha)
                alt=next(r for r in summaries if r['id']==cell['id'] and r['method']==METHODS[j] and r['scenario']=='moderate_binary_only' and r['alpha']==alpha)
                xx.append(null['calibrated_rate']);yy.append(alt['calibrated_rate'])
            ax.plot(xx,yy,marker='o',ms=3,color=color,label=SHORT[j])
        ax.set_xlabel('Held-out null rejection');ax.set_ylabel('Moderate binary-only power');ax.set_title(f"Varying dependence | n={cell['n']}");ax.legend(fontsize=9)
    fig.suptitle('SIMULATED data | operating curves from frozen calibration grid; no evaluation tuning')
    fig.tight_layout();save(fig,'power_vs_size.png')
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),sharey=True)
    for ax,scenario in zip(axes,['null','moderate_binary_only']):
        cell=cells()[1]
        for jj,j in enumerate(SHOW):
            r=lookup[cell['id'],scenario,j];b=bootlookup[cell['id'],scenario,j]
            ax.vlines(jj-.08,r['conditional_lower'],r['conditional_upper'],color='#333333',lw=2)
            ax.vlines(jj+.08,b['lower'],b['upper'],color='#b84b4b',lw=2)
            ax.scatter(jj,r['calibrated_rate'],color='black',s=12)
        ax.set_xticks(range(5),[SHORT[j] for j in SHOW],rotation=35,ha='right');ax.set_title(scenario.replace('_',' ').title());ax.set_ylabel('Calibrated rejection')
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([0],[0],color='#333333',label='Conditional Wilson'),Line2D([0],[0],color='#b84b4b',label='Calibration plus evaluation bootstrap')],loc='upper center',ncol=2,bbox_to_anchor=(.5,1.0))
    fig.suptitle('SIMULATED varying dependence | n=800',y=1.07);fig.tight_layout();save(fig,'calibration_uncertainty.png')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),sharey=True)
    for ax,scenario in zip(axes,['null','moderate_binary_only']):
        for j,color in zip(SHOW,COLORS):
            z=[r for r in banks if r['id']==1 and r['scenario']==scenario and r['method']==METHODS[j]]
            ax.plot([r['bank']+1 for r in z],[r['evaluation_rate'] for r in z],marker='o',color=color,label=SHORT[j])
        ax.set_xticks(range(1,6));ax.set_xlabel('Independent calibration subbank of 1,000');ax.set_ylabel('Held-out rejection');ax.set_title(scenario.replace('_',' ').title());ax.legend(fontsize=9)
    fig.suptitle('SIMULATED varying dependence | n=800 | shared held-out evaluations across subbanks');fig.tight_layout();save(fig,'subbank_sensitivity.png')
    prediction=read(ROOT/'results/information_predictions.csv')
    fig,axes=plt.subplots(2,2,figsize=(13,8),sharey=True)
    for ax,cell in zip(axes.flat,cells()):
        pp=[next(r for r in prediction if r['environment']==cell['environment'] and r['scenario']==s) for s in altnames]
        ax.plot(range(8),[float(r['joint_asymptotic_power']) for r in pp],color='#a94040',ls='--',label='Joint local prediction')
        ax.plot(range(8),[float(r['marginal_asymptotic_power']) for r in pp],color='#333333',ls=':',label='Marginal local prediction')
        for j,color in [(0,'#333333'),(3,'#377eaa'),(7,'#a94040')]:
            rr=[lookup[cell['id'],s,j] for s in altnames]
            ax.errorbar(np.arange(8)+(j/7-.5)*.12,[r['nominal_rate'] for r in rr],
                yerr=[[r['nominal_rate']-r['nominal_lower'] for r in rr],[r['nominal_upper']-r['nominal_rate'] for r in rr]],
                fmt='o',ms=3,capsize=1,color=color,label=SHORT[j]+' nominal')
        ax.set_xticks(range(8),labels,rotation=25);ax.set_ylabel('Nominal association power');ax.set_ylim(0,1.03);ax.set_title(f"n={cell['n']} | {cell['environment']}")
    axes[0,0].legend(fontsize=9);fig.suptitle('SIMULATED data | asymptotic information predictions and finite-sample nominal power')
    fig.tight_layout();save(fig,'information_predictions.png')
    import scipy
    result=dict(python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,platform=platform.platform(),
                execution='Native Windows; calibration started with three workers, reduced to two under memory pressure; source and design unchanged; completed chunks verified and resumed',bootstrap_replicates=cfg['bootstrap_replicates'])
    (ROOT/'results/software.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Seven figures and result tables generated',flush=True)
if __name__=='__main__':main()
