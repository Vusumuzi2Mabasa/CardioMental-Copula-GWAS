"""Count all secondary analyses, verify source/checkpoints and replay the polygenic model."""
from pathlib import Path
import sys,csv,json,gzip
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.append(str(ROOT/'v3/src'))
from prepare import dump,sha
from freeze import verify
from polygenic import simulate,tag,V
from polygenic_summary import load,cell
from reml import fit
from finemap_summary import load as load_fm
def main():
    lock=verify();assert json.loads((ROOT/'v3/results/audit.json').read_text())['passed']
    fmcal,fmeval=load_fm('cal'),load_fm('eval');assert len(fmcal)==2400 and len(fmeval)==3840
    ids=set()
    for r in fmcal+fmeval:
        key=tuple(r[k] for k in ['population','region','regime','phase','scenario','replicate','LD']);assert key not in ids;ids.add(key)
        assert 0<=r['max_pp4']<=1
    for path in (V/'results/finemap_chunks').glob('*.json'):
        if path.name=='software.json':continue
        m=json.loads(path.read_text());assert m['identity']==lock['identity']
        import hashlib
        assert hashlib.md5(path.with_suffix('.csv').read_bytes()).hexdigest()==m['output_md5']
    pip_count=0;sets=0
    ld_qc=[]
    from prepare import table
    for path in sorted((V/'work/finemap_inputs').glob('*.json')):
        d=json.loads(path.read_text());matrices=[np.array(d['Rref'])]
        for r in d['records']:matrices.extend([np.array(r['Rc']),np.array(r['Rb'])])
        minimum=min(float(np.linalg.eigvalsh(a).min()) for a in matrices)
        symmetry=max(float(np.max(abs(a-a.T))) for a in matrices);diagonal=max(float(np.max(abs(np.diag(a)-1))) for a in matrices)
        assert minimum>=.05-1e-8 and symmetry<1e-10 and diagonal<1e-10
        ld_qc.append(dict(group=path.stem,matrices=len(matrices),min_eigenvalue=minimum,max_asymmetry=symmetry,max_diagonal_error=diagonal))
    assert len(ld_qc)==156;table(V/'results/LD_validation.csv',ld_qc)
    for path in (V/'results/finemap_chunks').glob('*_pip.csv.gz'):
        with gzip.open(path,'rt') as f:
            for r in csv.DictReader(f):assert -1e-12<=float(r['pip'])<=1+1e-12;pip_count+=1
    for path in (V/'results/finemap_chunks').glob('*_sets.csv'):
        with path.open() as f:
            for r in csv.DictReader(f):
                markers=[int(x) for x in r['markers'].split(';')];assert len(markers)==int(r['size']) and all(0<=j<64 for j in markers);sets+=1
    cal,ev=load('cal'),load('eval');assert len(cal)==3200 and len(ev)==4800
    ids=set()
    for r in cal+ev:
        key=tuple(r[k] for k in ['population','n','residual_truth','genetic_truth','phase','replicate','mode']);assert key not in ids;ids.add(key)
        assert r['status']=='ok' or r['asymptotic_p']==1
    for path in (V/'results/polygenic_chunks').glob('*.json'):
        d=json.loads(path.read_text());assert d['identity']==lock['identity'] and d['sha256']==sha(path.with_suffix('.csv'))
    # Replays use an unselected interior replicate, separate from the R fixtures.
    errors=[]
    for pop in ['YRI','LWK']:
        for n in [600,1200]:
            with np.load(V/f'work/kernels/{pop}_n{n}.npz') as stored:panel={k:stored[k] for k in stored.files}
            for re in [0.,.4]:
                latent,observed,_=simulate(panel,pop,n,re,0.,'eval',37)
                for mode,y in [('liability',latent),('observed_binary',observed)]:
                    result,_,_=fit(panel['eigenvalues'],y)
                    old=next(r for r in ev if r['population']==pop and r['n']==n and r['residual_truth']==re and r['genetic_truth']==0 and r['replicate']==37 and r['mode']==mode)
                    error=max(abs(result[k]-old[k]) for k in ['genetic_correlation','genetic_covariance','LR']);assert error<1e-8;errors.append(error)
    rc=json.loads((V/'results/tests/R_reml_validation.json').read_text());assert rc['passed'] and rc['fixtures']>=33
    assert json.loads((V/'results/tests/R_finemap_validation.json').read_text())['passed']
    for file in ['finemap_audit.json','polygenic_audit.json']:assert json.loads((V/'results'/file).read_text())['passed']
    extension=json.loads((V/'protocol/validation_source_lock.json').read_text())
    import hashlib
    extid=hashlib.sha256((lock['identity']+sha(V/'protocol/validation_extension.json')+sha(V/'src/polygenic_validation.py')).encode()).hexdigest();assert extid==extension['identity']
    fc,fe=load('cal',True),load('eval',True);assert len(fc)==16000 and len(fe)==9600
    seen=set()
    for r in fc+fe:
        key=tuple(r[k] for k in ['population','n','residual_truth','genetic_truth','phase','replicate','mode']);assert key not in seen;seen.add(key)
        offset=1000000 if r['phase']=='cal' else 2000000
        assert r['seed_replicate']==r['replicate']+offset
    for path in (V/'results/polygenic_validation/polygenic_chunks').glob('*.json'):
        d=json.loads(path.read_text());assert d['identity']==extid and d['sha256']==sha(path.with_suffix('.csv'))
    assert json.loads((V/'results/polygenic_validation/polygenic_audit.json').read_text())['passed']
    with (V/'results/polygenic_validation/polygenic_rates.csv').open() as f:summ=list(csv.DictReader(f))
    cutoffs=json.loads((V/'results/polygenic_validation/polygenic_thresholds.json').read_text())['thresholds']
    for s in summ:
        x=[r for r in fe if r['population']==s['population'] and r['n']==int(s['n']) and r['residual_truth']==float(s['residual_truth']) and r['genetic_truth']==float(s['genetic_truth']) and r['mode']==s['mode']]
        assert len(x)==200
        good=[r for r in x if r['status']=='ok'];truth=float(s['genetic_truth'])
        assert abs(np.mean([r['genetic_correlation'] for r in good])-float(s['mean_rg']))<1e-14
        assert abs(np.sqrt(np.mean([(r['genetic_correlation']-truth)**2 for r in good]))-float(s['rmse_rg']))<1e-14
        key=tag(s['population'],int(s['n']),float(s['residual_truth']))+'_'+s['mode'];cut=cutoffs[key]
        calls=sum(cut['active'] and r['status']=='ok' and r['LR']>cut['LR_cutoff'] for r in x)/200
        assert calls==float(s['calibrated_rate'])
    # Replay four paired bootstrap rounds using an independent quantile implementation.
    with (V/'results/polygenic_validation/polygenic_bootstrap.csv').open() as f:bs=list(csv.DictReader(f))
    bmap={(int(r['bootstrap']),int(r['n']),float(r['genetic_truth']),r['mode']):float(r['rate']) for r in bs}
    cb={k:sorted([r for r in fc if cell(r)==k],key=lambda r:r['replicate']) for k in cutoffs}
    eb={(k,g):sorted([r for r in fe if cell(r)==k and r['genetic_truth']==g],key=lambda r:r['replicate']) for k in cutoffs for g in [-.5,0.,.5]}
    rng=np.random.default_rng(29103113)
    for bi in range(4):
        for n in [600,1200]:
            for g in [-.5,0.,.5]:
                counts={'liability':[],'observed_binary':[]}
                for pop in ['YRI','LWK']:
                    for re in [0.,.4]:
                        ii=rng.integers(1000,size=1000);jj=rng.integers(200,size=200)
                        for mode in counts:
                            key=tag(pop,n,re)+'_'+mode
                            cutoff=np.quantile([cb[key][i]['LR'] for i in ii],.95,method='higher')
                            counts[mode].append(np.mean([cutoffs[key]['active'] and eb[key,g][j]['status']=='ok' and eb[key,g][j]['LR']>cutoff for j in jj]))
                for mode,values in counts.items():assert abs(np.mean(values)-bmap[bi,n,g,mode])<1e-14
    # Independently replay one fresh validation cohort per population/size/residual stratum.
    fresh_errors=[]
    for pop in ['YRI','LWK']:
        for n in [600,1200]:
            with np.load(V/f'work/kernels/{pop}_n{n}.npz') as stored:panel={k:stored[k] for k in stored.files}
            for re in [0.,.4]:
                latent,observed,_=simulate(panel,pop,n,re,0.,'eval',2000037)
                for mode,y in [('liability',latent),('observed_binary',observed)]:
                    result,_,_=fit(panel['eigenvalues'],y)
                    old=next(r for r in fe if r['population']==pop and r['n']==n and r['residual_truth']==re and r['genetic_truth']==0 and r['replicate']==37 and r['mode']==mode)
                    error=max(abs(result[k]-old[k]) for k in ['genetic_correlation','genetic_covariance','LR']);assert error<1e-8;fresh_errors.append(error)
    result=dict(passed=True,source_identity=lock['identity'],protocol_commit='85b0bc4',fresh_validation_commit='dd84837',fresh_validation_identity=extid,finemap_analyses=6240,finemap_cohorts_reused_from_v3=3120,retained_variant_PIP_rows=pip_count,retained_credible_sets=sets,initial_polygenic_cohorts=4000,fresh_polygenic_cohorts=12800,polygenic_outcome_fit_records=33600,replayed_polygenic_fits=len(errors)+len(fresh_errors),max_replay_error=max(errors+fresh_errors),fresh_summary_rows_checked=len(summ),fresh_bootstrap_rounds_replayed=4,independent_R_fixtures=rc['fixtures'],scope='Synthetic analyses of association, causal-variant sharing and polygenic genetic correlation; no real shared genetic causes established')
    dump(V/'results/audit.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
