"""Matched-subset comparisons and integrated research figures for Levels 2 and 3."""
from pathlib import Path
import sys,csv,json,gzip
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];sys.path.append(str(ROOT/'v3/src'));V=ROOT/'v4'
from prepare import table
from analyse import load as load_v3,LABELS,wilson
from model import SCENARIOS,NO_SHARED
from finemap_summary import load as load_fm,tag as fm_tag
from polygenic_summary import load as load_poly
def read(path):
    with path.open() as f:return list(csv.DictReader(f))
def save(fig,name):
    (V/'figures').mkdir(exist_ok=True);fig.savefig(V/f'figures/{name}.png',dpi=170,bbox_inches='tight');plt.close(fig)
def main():
    plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False})
    banks=load_v3();fm=load_fm('eval');agg=[]
    for regime in ['ideal','realistic']:
        for ld in ['matched_reference','in_sample']:
            for s in SCENARIOS:
                numerators=dict(abf_conventional=0,abf_calibrated=0,susie_conventional=0,susie_calibrated=0)
                for pop in ['YRI','LWK']:
                    for region in range(3):
                        tag=f'{pop}_r{region}_{regime}'
                        cut=max(r['split_pp4'] for sc in NO_SHARED for r in banks[tag,'cal',sc][:20])
                        rows=banks[tag,'eval',s][:20];assert len(rows)==20
                        numerators['abf_conventional']+=sum(r['split_pp4']>.8 for r in rows)
                        numerators['abf_calibrated']+=sum(r['split_pp4']>cut for r in rows)
                        x=[r for r in fm if r['population']==pop and r['region']==region and r['regime']==regime and r['scenario']==s and r['LD']==ld];assert len(x)==20
                        ft=json.loads((V/'results/finemap_thresholds.json').read_text())['thresholds'][tag+'_'+ld]
                        numerators['susie_conventional']+=sum(r['status']=='ok' and r['max_pp4']>.8 for r in x)
                        numerators['susie_calibrated']+=sum(ft['active'] and r['status']=='ok' and r['max_pp4']>ft['cutoff'] for r in x)
                row=dict(regime=regime,LD=ld,scenario=s,attempted=120)
                for k,count in numerators.items():
                    lo,hi=wilson(count,120);row.update({k:count/120,k+'_conditional_lower':lo,k+'_conditional_upper':hi})
                agg.append(row)
    table(V/'results/finemap_matched_comparison.csv',agg)
    methods=['abf_conventional','susie_conventional','abf_calibrated','susie_calibrated'];colors=['#176B87','#B3541E','#3D7C47','#8A4F7D'];labels=['Single causal PP4 > 0.8','SuSiE any pair PP4 > 0.8','Single causal calibrated','SuSiE calibrated']
    fig,axs=plt.subplots(2,1,figsize=(7.4,8),constrained_layout=True)
    chosen=['null','cardio_only','mental_only','shared_same','shared_opposing','distinct_low','distinct_high','shared_plus_distinct']
    for ax,regime in zip(axs,['ideal','realistic']):
        for j,k in enumerate(methods):
            rows=[next(r for r in agg if r['regime']==regime and r['LD']=='matched_reference' and r['scenario']==s) for s in chosen]
            ax.bar(np.arange(8)+(j-1.5)*.18,[r[k] for r in rows],.18,color=colors[j],label=labels[j])
        ax.set_xticks(range(8),['No\neffects','Cardio\nonly','Mental\nonly','Shared\n+','Shared\n−','Distinct\nlow LD','Distinct\nhigh LD','Shared\nplus extra'],fontsize=10);ax.set_ylim(0,1.05);ax.set_title(regime.capitalize()+' phenotypes');ax.set_ylabel('Sharing declaration proportion')
    fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,1.10),ncol=2,fontsize=10);save(fig,'fine_mapping_comparison')
    example=V/'results/finemap_chunks/eval_YRI_r0_realistic_shared_plus_distinct_pip.csv.gz'
    with gzip.open(example,'rt') as f:pips=list(csv.DictReader(f))
    variants=[r for r in read(ROOT/'v3/reference/variants.csv') if r['region']=='0'];pos=np.array([int(r['position'])/1e6 for r in variants])
    fig,axs=plt.subplots(2,2,figsize=(7.4,6.5),constrained_layout=True)
    for i,ld in enumerate(['matched_reference','in_sample']):
        for j,trait in enumerate(['c','b']):
            a=sorted([r for r in pips if r['LD']==ld and r['trait']==trait and r['replicate']=='0'],key=lambda r:int(r['marker']));assert len(a)==64
            ax=axs[i,j];p=np.array([float(r['pip']) for r in a]);truth=np.array([r['true_causal']=='1' for r in a])
            ax.vlines(pos,0,p,color=colors[j],lw=1);ax.scatter(pos,p,color=colors[j],s=12)
            ax.scatter(pos[truth],np.full(truth.sum(),1.035),marker='*',color='#A16E00',s=90,label='Generating causal variant')
            ax.set_ylim(-.02,1.10);ax.set_ylabel('Inclusion probability');ax.set_xlabel('GRCh37 position Mb',fontsize=11);ax.set_xticks([20,20.05,20.1]);ax.ticklabel_format(axis='x',useOffset=False);ax.set_title(('Cardiometabolic' if trait=='c' else 'Mental health')+'\n'+ld.replace('_',' '),fontsize=11)
    fig.legend(*axs[0,0].get_legend_handles_labels(),fontsize=11,loc='upper center',bbox_to_anchor=(.5,1.10));save(fig,'fine_mapping_pip')
    base=V/'results/polygenic_validation';rates=read(base/'polygenic_rates.csv');pa=read(base/'polygenic_aggregate.csv');raw=load_poly('eval',True)
    fig,axs=plt.subplots(1,2,figsize=(7.4,4.5),constrained_layout=True)
    for ax,n in zip(axs,[600,1200]):
        for j,mode in enumerate(['liability','observed_binary']):
            values=[[r['genetic_correlation'] for r in raw if r['n']==n and r['mode']==mode and r['genetic_truth']==g and r['status']=='ok'] for g in [-.5,0.,.5]]
            bp=ax.boxplot(values,positions=np.arange(3)+(j-.5)*.25,widths=.20,patch_artist=True,showfliers=False)
            for box in bp['boxes']:box.set(facecolor=colors[j],alpha=.65)
            ax.plot([],[],color=colors[j],lw=6,label='Liability benchmark' if mode=='liability' else 'Observed binary')
        ax.scatter(range(3),[-.5,0,.5],marker='x',color='black',s=70,zorder=4,label='Generating value');ax.set_xticks(range(3),['−0.5','0','0.5']);ax.set_ylim(-1.03,1.03);ax.set_xlabel('Generating genetic correlation',fontsize=11);ax.set_ylabel('Estimated genetic correlation');ax.set_title(f'n = {n}   800 cohorts per value',fontsize=11)
    fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,1.13),ncol=3,fontsize=10);save(fig,'genetic_correlation')
    fig,axs=plt.subplots(1,2,figsize=(7.4,4.2),constrained_layout=True)
    for ax,n in zip(axs,[600,1200]):
        for j,(key,label) in enumerate([('mean_phenotypic_correlation','Observed phenotypic r'),('mean_rg','Estimated genetic r')]):
            values=[];ses=[]
            for re in [0.,.4]:
                a=[float(r[key]) for r in rates if int(r['n'])==n and float(r['genetic_truth'])==0 and float(r['residual_truth'])==re and r['mode']=='observed_binary'];values.append(np.mean(a))
                field='genetic_correlation' if key=='mean_rg' else 'phenotypic_correlation'
                groups=[[r[field] for r in raw if r['population']==pop and r['n']==n and r['residual_truth']==re and r['genetic_truth']==0 and r['mode']=='observed_binary' and r['status']=='ok'] for pop in ['YRI','LWK']]
                ses.append(np.sqrt(sum(np.var(g,ddof=1)/len(g) for g in groups)/4))
            ax.errorbar([0,.4],values,yerr=1.95996398454*np.array(ses),marker='o',color=colors[j],label=label,capsize=3)
        ax.axhline(0,color='gray',lw=.8);ax.set_xticks([0,.4]);ax.set_xlabel('Residual liability correlation',fontsize=11);ax.set_ylabel('Correlation');ax.set_title(f'Genetic truth = 0   n = {n}',fontsize=11);ax.set_ylim(-.06,.22)
    fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,1.13),ncol=2,fontsize=10);save(fig,'environmental_separation')
    fig,axs=plt.subplots(1,2,figsize=(7.4,4.2),constrained_layout=True)
    for ax,n in zip(axs,[600,1200]):
        for j,method in enumerate(['nominal_rate','calibrated_rate']):
            a=[next(r for r in pa if int(r['n'])==n and r['mode']==mode and float(r['genetic_truth'])==0) for mode in ['liability','observed_binary']]
            vals=[float(r[method]) for r in a];ax.bar(np.arange(2)+(j-.5)*.3,vals,.3,color=colors[j],label='Chi square p < 0.05' if j==0 else 'Independent calibration')
            if j==1:
                lo=np.array([float(r['bootstrap_lower']) for r in a]);hi=np.array([float(r['bootstrap_upper']) for r in a]);ax.errorbar(np.arange(2)+.15,vals,yerr=np.vstack([np.maximum(0,np.array(vals)-lo),np.maximum(0,hi-np.array(vals))]),fmt='none',ecolor='black',capsize=3)
        ax.axhline(.05,color='gray',ls='--',lw=1);ax.set_xticks([0,1],['Liability','Observed binary'],fontsize=11);ax.set_ylim(0,.12);ax.set_ylabel('False positive rate');ax.set_title(f'Fresh validation   n = {n}',fontsize=11)
    fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,1.13),ncol=2,fontsize=10);save(fig,'genetic_correlation_calibration')
    print('Level 2 and 3 figures generated')
if __name__=='__main__':main()
