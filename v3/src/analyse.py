"""Prespecified rate summaries, calibration-aware uncertainty and research figures."""
import csv,json,math
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import SCENARIOS,NO_SHARED,generate
from prepare import V,ROOT,dump,table
from run import cellid,identity

METHODS=['Joint association','Tagging conjunction','Tagging calibrated','Coloc PP4 > 0.8','Coloc calibrated','Overlap PP4 > 0.8','Prior low PP4 > 0.8','Prior high PP4 > 0.8']
LABELS={'null':'No genetic effect','cardio_only':'Cardiometabolic only','mental_only':'Mental health only','shared_same':'Shared concordant','shared_opposing':'Shared opposing','distinct_low':'Distinct low LD','distinct_high':'Distinct high LD','shared_plus_distinct':'Shared plus distinct'}
def load():
    banks=defaultdict(list)
    for path in sorted((V/'results/chunks').rglob('*.csv')):
        with path.open() as f:
            for row in csv.DictReader(f):
                for k,v in list(row.items()):
                    if k not in ['population','regime','phase','scenario','fingerprint','status','copula_status']:row[k]=float(v)
                banks[(cellid(row['population'],int(row['region']),row['regime']),row['phase'],row['scenario'])].append(row)
    for key in banks:banks[key].sort(key=lambda r:r['replicate'])
    return banks
def arr(rows,key):return np.array([r[key] for r in rows])
def thresholds_from(banks,tag,indices=None):
    n=len(banks[tag,'cal','null']);k=int(np.floor(.05*(n+1)))
    a={s:banks[tag,'cal',s] if indices is None else [banks[tag,'cal',s][i] for i in indices[s]] for s in NO_SHARED}
    return dict(split_pp4=max(np.sort(arr(a[s],'split_pp4'))[n-k] for s in NO_SHARED),iut_p=min(np.sort(arr(a[s],'iut_p'))[k-1] for s in NO_SHARED[:3]))
def calls(rows,t):
    return np.column_stack([arr(rows,'joint_p')<.05,arr(rows,'iut_p')<.05,arr(rows,'iut_p')<t['iut_p'],arr(rows,'split_pp4')>.8,arr(rows,'split_pp4')>t['split_pp4'],arr(rows,'overlap_pp4')>.8,arr(rows,'prior_low_pp4')>.8,arr(rows,'prior_high_pp4')>.8])
def wilson(k,n):
    z=1.95996398454;p=k/n;center=(p+z*z/(2*n))/(1+z*z/n);half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return float(center-half),float(center+half)
def main():
    banks=load();saved=json.loads((V/'results/thresholds.json').read_text());assert saved['identity']==identity();ts=saved['thresholds']
    cfg=json.loads((V/'protocol/design.json').read_text());summ=[];agg=defaultdict(list);phen=[];cop=[]
    for tag,t in ts.items():
        pop,reg,regime=tag.split('_');region=int(reg[1:])
        recompute=thresholds_from(banks,tag)
        assert all(abs(t[k]-recompute[k])<1e-15 for k in recompute)
        for sc in SCENARIOS:
            rows=banks[tag,'eval',sc];assert len(rows)==cfg['evaluation_replicates'];cc=calls(rows,t)
            for j,method in enumerate(METHODS):
                k=int(cc[:,j].sum());lo,hi=wilson(k,len(rows))
                summ.append(dict(population=pop,region=region,regime=regime,scenario=sc,method=method,attempted=len(rows),discoveries=k,rate=k/len(rows),lower=lo,upper=hi))
            agg[regime,sc].append(cc)
            phen.append(dict(population=pop,region=region,regime=regime,scenario=sc,**{k:float(np.mean(arr(rows,k))) for k in ['n_complete','prevalence','true_prevalence','pressure_mean','pressure_sd','treated_fraction','missing_fraction','comorbidity','comorbidity_or']}))
            chosen=[r for r in rows if r['copula_status']!='not_planned'];n=len(chosen)
            cop.append(dict(population=pop,region=region,regime=regime,scenario=sc,attempted=n,failures=sum(r['copula_status']!='ok' for r in chosen),copula_nominal_rate=float(np.mean(arr(chosen,'copula_p')<.05)),marginal_anchor_joint_rate=float(np.mean(arr(chosen,'anchor_joint')<.05)),tagging_anchor_rate=float(np.mean(arr(chosen,'anchor_iut')<.05))))
    table(V/'results/rates.csv',summ);table(V/'results/phenotypes.csv',phen);table(V/'results/copula_diagnostic.csv',cop)
    # Stratified paired resampling includes uncertainty from composite calibration banks.
    rng=np.random.default_rng(29103103);B=400;draws=np.zeros((B,2,8,len(METHODS)))
    for bi in range(B):
        for ri,regime in enumerate(['ideal','realistic']):
            samples={s:[] for s in SCENARIOS}
            for tag in [k for k in ts if k.endswith('_'+regime)]:
                n=cfg['calibration_replicates'];idx={s:rng.integers(n,size=n) for s in NO_SHARED}
                t=thresholds_from(banks,tag,idx)
                for s in SCENARIOS:
                    n=cfg['evaluation_replicates'];ii=rng.integers(n,size=n);rr=banks[tag,'eval',s]
                    samples[s].append(calls([rr[i] for i in ii],t))
            for si,s in enumerate(SCENARIOS):draws[bi,ri,si]=np.vstack(samples[s]).mean(0)
        if (bi+1)%100==0:print('bootstrap',bi+1,flush=True)
    np.savez_compressed(V/'results/bootstrap_draws.npz',rates=draws)
    aggregate=[];contrasts=[]
    for ri,regime in enumerate(['ideal','realistic']):
        for si,s in enumerate(SCENARIOS):
            rates=np.vstack(agg[regime,s]).mean(0)
            for j,method in enumerate(METHODS):
                lo,hi=np.quantile(draws[:,ri,si,j],[.025,.975])
                aggregate.append(dict(regime=regime,scenario=s,method=method,attempted=1200,rate=float(rates[j]),bootstrap_lower=float(lo),bootstrap_upper=float(hi)))
            for j,k in [(4,3),(5,3),(1,3)]:
                lo,hi=np.quantile(draws[:,ri,si,j]-draws[:,ri,si,k],[.025,.975])
                contrasts.append(dict(regime=regime,scenario=s,contrast=METHODS[j]+' minus '+METHODS[k],difference=float(rates[j]-rates[k]),bootstrap_lower=float(lo),bootstrap_upper=float(hi)))
    table(V/'results/aggregate_rates.csv',aggregate);table(V/'results/paired_contrasts.csv',contrasts)
    truth=[]
    for pop in ['YRI','LWK']:
        for region in range(3):
            h=np.load(V/f'reference/region_{region}.npz')[pop+'_study'].astype(float)
            meta=json.loads((V/f'reference/region_{region}.json').read_text())[pop]
            sigma=2*np.cov(h,rowvar=False,ddof=0)
            for sc in SCENARIOS:
                *_,bc,bm=generate(np.random.default_rng(991),h,meta,sc,'ideal',n=8)
                total=float(bc@sigma@bm);diagonal=float(np.sum(bc*bm*np.diag(sigma)))
                truth.append(dict(population=pop,region=region,scenario=sc,shared_direct_variants=int(np.sum((bc!=0)&(bm!=0))),genetic_predictor_covariance=total,same_variant_diagonal=diagonal,between_variant_LD_component=total-diagonal))
    table(V/'results/generating_genetic_covariance.csv',truth)
    figures(summ,aggregate,phen,banks,ts)

def save(fig,name):
    (V/'figures').mkdir(exist_ok=True);fig.savefig(V/f'figures/{name}.png',dpi=170,bbox_inches='tight');plt.close(fig)
def figures(summ,agg,phen,banks,ts):
    plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False})
    short=['No\neffects','Cardio\nonly','Mental\nonly','Shared\n+','Shared\n−','Distinct\nlow LD','Distinct\nhigh LD','Shared\nplus extra']
    colors=['#176B87','#B3541E','#3D7C47','#8A4F7D'];pair=list(csv.DictReader((V/'reference/causal_pairs.csv').open()))
    fig,axs=plt.subplots(3,2,figsize=(7.4,8.1),constrained_layout=True)
    for row,ax in zip(pair,axs.flat):
        region=int(row['region']);pop=row['population'];r=np.load(V/f'reference/region_{region}.npz')[pop+'_study_R']**2
        im=ax.imshow(r,vmin=0,vmax=1,cmap='viridis',origin='lower');ax.set_title(f'{pop} region {region+1}  pair r²={float(row["study_r2_high"]):.2f}',fontsize=11);ax.set_xlabel('Marker index');ax.set_ylabel('Marker index')
    fig.colorbar(im,ax=axs.ravel().tolist(),label='Haplotype r²',shrink=.65)
    fig.suptitle('African haplotype LD   GRCh37 chromosome 22',fontsize=13);save(fig,'african_ld')
    for regime in ['ideal','realistic']:
        fig,ax=plt.subplots(figsize=(7.4,4.8),constrained_layout=True);x=np.arange(8);width=.19
        for j,method in enumerate([METHODS[0],METHODS[1],METHODS[3],METHODS[4]]):
            rows=[next(r for r in agg if r['regime']==regime and r['scenario']==s and r['method']==method) for s in SCENARIOS]
            y=np.array([r['rate'] for r in rows]);lo=np.array([r['bootstrap_lower'] for r in rows]);hi=np.array([r['bootstrap_upper'] for r in rows])
            ax.bar(x+(j-1.5)*width,y,width,color=colors[j],label=method)
            ax.errorbar(x+(j-1.5)*width,y,yerr=np.vstack([np.maximum(0,y-lo),np.maximum(0,hi-y)]),fmt='none',ecolor='#333333',capsize=2,lw=.8)
        ax.set_xticks(x,short,fontsize=10);ax.set_ylim(0,1.05);ax.set_ylabel('Discovery proportion');ax.set_title(regime.capitalize()+' phenotypes   1200 datasets per architecture',fontsize=12);ax.legend(ncol=2,loc='upper center',bbox_to_anchor=(.5,1.31),fontsize=10);save(fig,'discovery_'+regime)
    fig,axs=plt.subplots(1,2,figsize=(7.4,4.4),constrained_layout=True)
    for ax,regime in zip(axs,['ideal','realistic']):
        for j,method in enumerate([METHODS[3],METHODS[4]]):
            for sc,mark in [('shared_same','o'),('distinct_high','s')]:
                rows=[r for r in summ if r['regime']==regime and r['scenario']==sc and r['method']==method]
                x=np.arange(6)+(.10 if j else -.10)
                ax.plot(x,[r['rate'] for r in rows],marker=mark,label=('Calibrated' if j else 'Fixed cutoff')+' / '+('shared' if sc=='shared_same' else 'distinct'),ls='-' if j else '--',color=colors[j+2],alpha=.8)
        ax.axhline(.05,color='gray',lw=1);ax.set_xticks(range(6),['YRI1','YRI2','YRI3','LWK1','LWK2','LWK3']);ax.set_ylim(-.01,1);ax.set_title(regime.capitalize());ax.set_ylabel('Shared causal variant declaration')
    fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,1.17),ncol=2,fontsize=11);save(fig,'identifiability')
    fig,axs=plt.subplots(1,2,figsize=(7.4,4.4),constrained_layout=True)
    for ax,regime in zip(axs,['ideal','realistic']):
        for j,(method,label) in enumerate([(METHODS[6],'p12 = 10⁻⁶'),(METHODS[3],'p12 = 10⁻⁵'),(METHODS[7],'p12 = 10⁻⁴'),(METHODS[5],'Overlapping subjects')]):
            ss=['shared_same','shared_opposing','distinct_low','distinct_high'];rr=[next(r for r in agg if r['regime']==regime and r['scenario']==s and r['method']==method) for s in ss]
            ax.plot(range(4),[r['rate'] for r in rr],marker='o',color=colors[j],label=label)
        ax.set_xticks(range(4),['Shared +','Shared −','Distinct low','Distinct high'],rotation=20);ax.set_title(regime.capitalize());ax.set_ylim(0,1);ax.set_ylabel('PP4 > 0.8 proportion')
    fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,1.17),ncol=2,fontsize=10);save(fig,'prior_overlap_sensitivity')
    fig,axs=plt.subplots(1,3,figsize=(7.4,3.4),constrained_layout=True)
    for ax,key,title in zip(axs,['prevalence','missing_fraction','comorbidity_or'],['Observed depression\nprevalence','Missing observation\nfraction','Phenotypic comorbidity\nodds ratio']):
        for j,regime in enumerate(['ideal','realistic']):
            a=[r[key] for r in phen if r['regime']==regime];ax.boxplot(a,positions=[j],widths=.5,patch_artist=True,boxprops=dict(facecolor=colors[j],alpha=.6))
        ax.set_xticks([0,1],['Ideal','Realistic'],fontsize=11);ax.set_title(title,fontsize=10.5)
    fig.suptitle('Synthetic phenotype characteristics across design cells',fontsize=12);save(fig,'phenotype_checks')
    fig,axs=plt.subplots(1,2,figsize=(7.4,4.4),constrained_layout=True)
    for ax,regime in zip(axs,['ideal','realistic']):
        for j,s in enumerate(NO_SHARED):
            vals=np.concatenate([arr(banks[tag,'eval',s],'split_pp4') for tag in ts if tag.endswith('_'+regime)])
            ax.ecdf(vals,label={'null':'No effects','cardio_only':'Cardio only','mental_only':'Mental only','distinct_low':'Distinct low LD','distinct_high':'Distinct high LD'}[s])
        ax.axvline(.8,color='gray',ls='--');ax.set_xlabel('Split-sample PP4');ax.set_ylabel('Held-out cumulative proportion');ax.set_title(regime.capitalize())
    fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,1.20),ncol=3,fontsize=11);save(fig,'composite_null')
if __name__=='__main__':main()
