"""Summaries and scientific figures from audited simulation records."""
from pathlib import Path
import json,csv,itertools
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from benchmark import cells,ROOT,write_csv
from model import METHODS

SHORT=dict(zip(METHODS,['Bonf','MargCov','Indep','C-model','C-robust','V-model','V-robust']))
COLORS=dict(zip(METHODS,['#67737c','#202020','#d97732','#3975a6','#3b9d96','#9a568f','#b63c42']))
SCENARIOS=['continuous_only','binary_only','shared','opposing']
LABELS=['Continuous only','Binary only','Shared','Opposing']


def read(path):
    with Path(path).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
def mean_se(a):
    a=np.asarray(a,dtype=float)
    return (float(np.mean(a)),float(np.std(a,ddof=1)/np.sqrt(len(a)))) if len(a)>1 else (float(np.mean(a)),float('nan'))
def binomial_ci(p,n):
    z=1.95996398454
    den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return max(0.,center-half),min(1.,center+half)
def errors(rows,field):
    return np.maximum(0,np.array([[r[field]-r[field+'_lower'] for r in rows],
                                  [r[field+'_upper']-r[field] for r in rows]]))
def save(fig,name):
    fig.savefig(ROOT/'figures'/name,dpi=220,bbox_inches='tight',facecolor='white');plt.close(fig)
def metadata(c):return {k:c[k] for k in ['id','group','setting','n','environment','prevalence','maf','ld']}


def main():
    assert json.loads((ROOT/'results/audit.json').read_text())['passed']
    plt.rcParams.update({'font.size':14,'axes.spines.top':False,'axes.spines.right':False})
    summaries=[];effects=[];pairs=[];null_curves=[];all_rows=[]
    for c in cells():
        folder=ROOT/'results/cells'/f"cell_{c['id']:02d}"
        rows=read(folder/'scores.csv');fits=read(folder/'null_fits.csv');ee=read(folder/'effects.csv')
        freeze=json.loads((folder/'threshold_freeze.json').read_text())
        all_rows.extend((c,r) for r in rows)
        for scenario,method in itertools.product(['null']+SCENARIOS,METHODS):
            selected=[r for r in rows if r['scenario']==scenario and r['method']==method]
            s=dict(**metadata(c),scenario=scenario,method=method,attempted=len(selected),
                valid=sum(r['status']=='ok' for r in selected),failures=sum(r['status']!='ok' for r in selected),
                calibration_failures=freeze['calibration_failures'][method],threshold=freeze['thresholds'][method],
                mean_seconds=float(np.mean([float(r['seconds']) for r in selected])))
            for field in ['nominal_05','nominal_01','calibrated_05','panel_FWER_05']:
                vals=[float(r[field]) for r in selected if r[field]!='']
                if vals and (field!='panel_FWER_05' or scenario=='null'):
                    s[field],s[field+'_mcse']=mean_se(vals)
                else:s[field]=s[field+'_mcse']=float('nan')
                if vals and np.isfinite(s[field]):
                    s[field+'_lower'],s[field+'_upper']=binomial_ci(s[field],len(vals))
                else:s[field+'_lower']=s[field+'_upper']=float('nan')
            summaries.append(s)
        comparisons=[('Varying sandwich score','Constant sandwich score'),
                     ('Varying sandwich score','Marginal covariance-aware'),
                     ('Varying model score','Constant model score'),
                     ('Varying model score','Marginal Bonferroni')]
        lookup={(r['scenario'],r['method'],int(r['replicate'])):r for r in rows}
        for scenario,(a,b),metric in itertools.product(['null']+SCENARIOS,comparisons,['nominal_05','calibrated_05']):
            n=c['null'] if scenario=='null' else c['alt'];diffs=[]
            for rep in range(n):
                aa=lookup[scenario,a,rep][metric];bb=lookup[scenario,b,rep][metric]
                if aa!='' and bb!='':diffs.append(float(aa)-float(bb))
            m,se=mean_se(diffs) if diffs else (float('nan'),float('nan'))
            pairs.append(dict(**metadata(c),scenario=scenario,method_a=a,method_b=b,metric=metric,
                attempted_pairs=n,valid_pairs=len(diffs),difference=m,mcse=se,lower=m-1.96*se,upper=m+1.96*se))
        for scenario,model,trait in itertools.product(['null']+SCENARIOS,['Constant','Varying'],['continuous','binary_liability']):
            z=[r for r in ee if r['scenario']==scenario and r['model']==model and r['trait']==trait]
            good=[r for r in z if r['status']=='ok']
            err=np.array([float(r['estimate'])-float(r['truth']) for r in good])
            out=dict(**metadata(c),scenario=scenario,model=model,trait=trait,attempted=len(z),valid=len(good),
                bias=float(np.mean(err)) if len(err) else float('nan'),rmse=float(np.sqrt(np.mean(err*err))) if len(err) else float('nan'))
            for kind in ['model','sandwich']:
                out[kind+'_coverage'],out[kind+'_coverage_mcse']=mean_se([int(r[kind+'_coverage']) for r in z])
                out[kind+'_valid_coverage']=float(np.mean([int(r[kind+'_coverage']) for r in good])) if good else float('nan')
                out[kind+'_mean_se']=float(np.mean([float(r[kind+'_se']) for r in good])) if good else float('nan')
                out[kind+'_coverage_lower'],out[kind+'_coverage_upper']=binomial_ci(out[kind+'_coverage'],len(z))
            effects.append(out)
        if c['group']=='primary' and c['environment']!='t_varying':
            for model in ['Constant','Varying']:
                good=[r for r in fits if r['scenario']=='null' and r['model']==model and r['status']=='ok']
                null_curves.append((c,model,good))
    write_csv(ROOT/'results/score_summary.csv',summaries)
    write_csv(ROOT/'results/paired_contrasts.csv',pairs)
    write_csv(ROOT/'results/effect_summary.csv',effects)
    sl={(r['id'],r['scenario'],r['method']):r for r in summaries}
    byid={c['id']:c for c in cells()}

    # Nominal versus empirically calibrated null rejection.
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharey=True)
    focus=['Marginal covariance-aware','Constant model score','Constant sandwich score','Varying model score','Varying sandwich score']
    envs=['independent','constant','varying','t_varying']
    for ax,env in zip(axes.flat,envs):
        for j,n in enumerate([800,1600]):
            c=next(c for c in cells() if c['group']=='primary' and c['environment']==env and c['n']==n)
            for k,m in enumerate(focus):
                r=sl[c['id'],'null',m];x=k+(j-.5)*.23
                ax.errorbar(x,r['nominal_05'],yerr=errors([r],'nominal_05'),fmt='o' if j==0 else 's',
                    color=COLORS[m],capsize=2,label=f'n={n}' if k==0 else None)
        ax.axhline(.05,color='gray',ls='--');ax.set_title(env.replace('_',' ').title())
        ax.set_xticks(range(len(focus)),[SHORT[m] for m in focus],rotation=25)
        ax.set_ylabel('Nominal SNP null rejection');ax.legend(fontsize=12)
    fig.suptitle('SIMULATED data | nominal alpha 0.05 | bars: 95% binomial Wilson intervals')
    fig.tight_layout();save(fig,'null_calibration.png')

    # Nominal and calibrated power on identical datasets at n=800.
    focus=['Marginal Bonferroni','Marginal covariance-aware','Constant sandwich score','Varying sandwich score']
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharey=True)
    for row,env in enumerate(['constant','varying']):
        c=next(c for c in cells() if c['group']=='primary' and c['environment']==env and c['n']==800)
        for col,metric in enumerate(['nominal_05','calibrated_05']):
            ax=axes[row,col]
            for m in focus:
                rs=[sl[c['id'],s,m] for s in SCENARIOS]
                ax.errorbar(range(4),[r[metric] for r in rs],yerr=errors(rs,metric),
                    marker='o',capsize=2,color=COLORS[m],label=SHORT[m])
            ax.set_title(f"{env.title()} dependence | {'nominal' if col==0 else 'calibrated'}")
            ax.set_xticks(range(4),LABELS,rotation=15);ax.set_ylim(0,1.05);ax.set_ylabel('Focal association power')
            ax.legend(fontsize=12,loc='lower right')
    fig.suptitle('SIMULATED data | n=800 | calibrated results condition on frozen null banks')
    fig.tight_layout();save(fig,'association_power.png')

    # Recovery of generating Gaussian copula dependence under the global null.
    grid=np.linspace(-2.5,2.5,100)
    fig,axes=plt.subplots(1,3,figsize=(13,4.2),sharey=True)
    for ax,env in zip(axes,['independent','constant','varying']):
        true=np.zeros_like(grid) if env=='independent' else np.full_like(grid,.4) if env=='constant' else np.tanh(np.arctanh(.3)+.6*grid)
        ax.plot(grid,true,color='black',ls='--',label='Generating curve')
        for c,model,good in null_curves:
            if c['environment']!=env or c['n']!=800:continue
            a=np.array([[float(r['rho_intercept']),float(r['rho_age'])] for r in good])
            curves=np.tanh(a[:,0,None]+a[:,1,None]*grid)
            avg=curves.mean(0);se=curves.std(0,ddof=1)/np.sqrt(len(curves))
            color=COLORS[model+' sandwich score']
            ax.plot(grid,avg,color=color,label=model)
            ax.fill_between(grid,avg-1.96*se,avg+1.96*se,color=color,alpha=.2)
        ax.set_title(env.title());ax.set_xlabel('Age-like covariate');ax.set_ylim(-1,1);ax.legend(fontsize=12)
    axes[0].set_ylabel('Gaussian copula correlation')
    fig.suptitle('SIMULATED Gaussian-copula null datasets | n=800 | mean fitted curves')
    fig.tight_layout();save(fig,'dependence_recovery.png')

    # Focal-SNP QQ plots avoid treating the eight markers as independent replicates.
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    for ax,env in zip(axes,['varying','t_varying']):
        c=next(c for c in cells() if c['group']=='primary' and c['environment']==env and c['n']==800)
        for m in ['Marginal covariance-aware','Constant model score','Varying model score','Varying sandwich score']:
            ps=np.sort([float(r['p']) for cc,r in all_rows if cc['id']==c['id'] and r['scenario']=='null' and r['method']==m and r['status']=='ok'])
            expected=-np.log10((np.arange(len(ps))+.5)/len(ps))
            ax.plot(expected,-np.log10(np.maximum(ps,1e-300)),color=COLORS[m],label=SHORT[m])
        ax.plot([0,3],[0,3],color='gray',ls='--');ax.set_title(env.replace('_',' ').title())
        ax.set_xlabel('Expected -log10(P)');ax.set_ylabel('Observed -log10(P)');ax.legend(fontsize=12)
    fig.suptitle('SIMULATED null data | focal SNP | diagnostic, not GWAS-tail validation')
    fig.tight_layout();save(fig,'score_qq.png')

    fig,axes=plt.subplots(1,2,figsize=(11,4.5),sharey=True)
    c=next(c for c in cells() if c['group']=='primary' and c['environment']=='varying' and c['n']==800)
    for ax,trait in zip(axes,['continuous','binary_liability']):
        for model in ['Constant','Varying']:
            for kind in ['model','sandwich']:
                z=[next(r for r in effects if r['id']==c['id'] and r['scenario']==s and r['model']==model and r['trait']==trait) for s in SCENARIOS]
                ax.errorbar(range(4),[r[kind+'_coverage'] for r in z],yerr=errors(z,kind+'_coverage'),
                    marker='o',capsize=2,ls='--' if kind=='model' else '-',color=COLORS[model+' sandwich score'],label=model+' '+kind)
        ax.axhline(.95,color='gray',ls=':');ax.set_xticks(range(4),LABELS,rotation=15)
        ax.set_title(trait.replace('_',' ').title());ax.set_ylim(0,1.1);ax.set_ylabel('Generating-effect interval coverage');ax.legend(fontsize=12)
    fig.suptitle('SIMULATED varying dependence | n=800 | 30 prespecified full fits per architecture')
    fig.tight_layout();save(fig,'effect_interval_coverage.png')

    fig,axes=plt.subplots(1,2,figsize=(12,4.5),sharey=True)
    settings=['primary','low_prevalence','high_prevalence','low_maf','LD']
    for ax,env in zip(axes,['constant','varying']):
        for m in ['Marginal covariance-aware','Constant sandwich score','Varying sandwich score']:
            z=[]
            for setting in settings:
                c=next(c for c in cells() if c['n']==800 and c['environment']==env and c['setting']==setting)
                z.append(sl[c['id'],'binary_only',m])
            ax.errorbar(range(len(settings)),[r['calibrated_05'] for r in z],
                yerr=errors(z,'calibrated_05'),marker='o',capsize=2,color=COLORS[m],label=SHORT[m])
        ax.set_xticks(range(len(settings)),['Main','Prev .05','Prev .30','MAF .05','LD'],rotation=20)
        ax.set_title(env.title()+' dependence');ax.set_ylabel('Calibrated binary-only signal power');ax.set_ylim(0,1.05);ax.legend(fontsize=12)
    fig.suptitle('SIMULATED sensitivity analysis | n=800 | independent calibration per setting')
    fig.tight_layout();save(fig,'sensitivity_power.png')

    fig,ax=plt.subplots(figsize=(8,4.5))
    runtime=[]
    for n in [800,1600]:
        for m in ['Marginal covariance-aware','Constant model score','Varying model score']:
            z=[r for r in summaries if r['group']=='primary' and r['n']==n and r['method']==m and r['scenario']=='null']
            runtime.append(dict(n=n,method=m,seconds=np.mean([r['mean_seconds'] for r in z])))
    for m in ['Marginal covariance-aware','Constant model score','Varying model score']:
        zz=[r for r in runtime if r['method']==m]
        ax.plot([r['n'] for r in zz],[r['seconds'] for r in zz],marker='o',label=SHORT[m],color=COLORS[m])
    ax.set_xlabel('People per dataset');ax.set_ylabel('Null fit + eight-SNP scan (seconds)')
    ax.set_title('SIMULATED workloads | two concurrent workers | hardware-specific timings');ax.legend()
    fig.tight_layout();save(fig,'runtime.png')
    write_csv(ROOT/'results/runtime_summary.csv',runtime)
    print('Summary tables and seven scientific figures generated')

if __name__=='__main__':main()
