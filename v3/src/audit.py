"""Independent record, preservation, replay and summary checks before publication."""
import argparse,csv,json,hashlib,subprocess,time
import numpy as np
from scipy.stats import binomtest
from prepare import V,ROOT,sha,dump
from model import SCENARIOS,NO_SHARED,generate,analyse_dataset
from run import identity,cellid
from analyse import load,calls,METHODS,thresholds_from

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--summaries',action='store_true');args=ap.parse_args()
    cfg=json.loads((V/'protocol/design.json').read_text());fp=identity();assert json.loads((V/'protocol/source_lock.json').read_text())['identity']==fp
    # Byte preservation includes all old scientific outputs, figures, engine and Word reports.
    excluded={'README.md','workflow.txt','checksums_public.sha256','.gitignore','.github/workflows/numerical-smoke.yml'}
    old=ROOT/'baseline/checksums_v2.sha256';preserved=0
    for line in old.read_text().splitlines():
        digest,path=line.split('  ',1)
        if path in excluded:continue
        assert sha(ROOT/path)==digest,'Version two changed: '+path
        preserved+=1
    assert sha(ROOT/'baseline/README_v2.md')==hashlib.sha256(subprocess.check_output(['git','-c','core.fsmonitor=false','show','a648ebd:README.md'],cwd=ROOT)).hexdigest()
    banks=load();n=0;unique=set();failures=0;copfail=0;plannedcop=0
    threshold_path=V/'results/thresholds.json';cut=json.loads(threshold_path.read_text());assert cut['identity']==fp
    previous=json.loads((V/'results/audit.json').read_text()) if (V/'results/audit.json').exists() else {}
    saved_timeline=previous.get('calibration_before_evaluation',False) and previous.get('source_identity')==fp and previous.get('threshold_sha256')==sha(threshold_path)
    for path in (V/'results/chunks').rglob('*.csv'):
        meta=json.loads(path.with_suffix('.json').read_text());assert meta['identity']==fp and meta['sha256']==sha(path)
        if path.name.startswith('eval_') and not saved_timeline:assert path.stat().st_mtime>=threshold_path.stat().st_mtime,'Evaluation precedes calibration freeze'
    for pop in cfg['populations']:
        for region in range(3):
            for regime in cfg['regimes']:
                tag=cellid(pop,region,regime)
                t=thresholds_from(banks,tag)
                assert all(t[k]==cut['thresholds'][tag][k] for k in t)
                for phase,scenarios in [('cal',NO_SHARED),('eval',SCENARIOS)]:
                    for scenario in scenarios:
                        rows=banks[tag,phase,scenario];assert len(rows)==200
                        assert [int(r['replicate']) for r in rows]==list(range(200))
                        for r in rows:
                            n+=1;assert r['fingerprint'] not in unique;unique.add(r['fingerprint'])
                            assert r['true_shared']==int(scenario.startswith('shared'))
                            assert 2700<=r['n_complete']<=3200
                            assert all(0<=r[k]<=1 for k in ['iut_p','joint_p','split_pp4','split_pp3','overlap_pp4','prior_low_pp4','prior_high_pp4'])
                            failures+=r['status']!='ok'
                            if r['copula_status']!='not_planned':plannedcop+=1;copfail+=r['copula_status']!='ok'
    assert n==31200 and len(unique)==n
    # Replay one independently selected evaluation architecture for every design cell.
    replay=[]
    for pi,pop in enumerate(cfg['populations']):
        for region in range(3):
            for ri,regime in enumerate(cfg['regimes']):
                si=(pi*6+region*2+ri)%8;sc=SCENARIOS[si];rep=37
                seed=[29103102,pi,region,ri,si,rep]
                h=np.load(V/f'reference/region_{region}.npz')[pop+'_study'];meta=json.loads((V/f'reference/region_{region}.json').read_text())[pop]
                X,y,d,G,_,_,_=generate(np.random.default_rng(np.random.SeedSequence(seed)),h,meta,sc,regime)
                got,_=analyse_dataset(X,y,d,G,meta,False);old=banks[cellid(pop,region,regime),'eval',sc][rep]
                err=max(abs(got[k]-old[k]) for k in got if k not in ['copula_status'])
                assert err<1e-12;replay.append(err)
    result=dict(passed=True,datasets=n,calibration_datasets=12000,evaluation_datasets=19200,unique_fingerprints=len(unique),main_failures=failures,copula_diagnostic_datasets=plannedcop,copula_failures=copfail,replayed_datasets=len(replay),max_replay_error=max(replay),preserved_old_files=preserved,calibration_before_evaluation=True,threshold_sha256=sha(threshold_path),timeline_evidence='Initial execution file ordering, retained by threshold hash for portable replay' if saved_timeline else 'Initial execution file ordering checked before publication',source_identity=fp,protocol_commit='5999591',audited_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    if args.summaries:
        with (V/'results/rates.csv').open() as f:rates=list(csv.DictReader(f))
        for r in rates:
            rows=banks[cellid(r['population'],int(r['region']),r['regime']),'eval',r['scenario']];t=cut['thresholds'][cellid(r['population'],int(r['region']),r['regime'])]
            j=METHODS.index(r['method']);k=int(calls(rows,t)[:,j].sum());n=len(rows)
            assert k==int(r['discoveries']) and k/n==float(r['rate'])
            ci=binomtest(k,n).proportion_ci(method='wilson')
            assert abs(ci.low-float(r['lower']))<1e-12 and abs(ci.high-float(r['upper']))<1e-12
        with (V/'results/aggregate_rates.csv').open() as f:agg=list(csv.DictReader(f))
        d=np.load(V/'results/bootstrap_draws.npz')['rates'];assert d.shape==(400,2,8,8)
        for r in agg:
            sub=[x for x in rates if all(x[k]==r[k] for k in ['regime','scenario','method'])]
            assert len(sub)==6 and abs(sum(int(x['discoveries']) for x in sub)/1200-float(r['rate']))<1e-15
            ri=cfg['regimes'].index(r['regime']);si=SCENARIOS.index(r['scenario']);j=METHODS.index(r['method'])
            lo,hi=np.quantile(d[:,ri,si,j],[.025,.975]);assert abs(lo-float(r['bootstrap_lower']))<1e-15 and abs(hi-float(r['bootstrap_upper']))<1e-15
        result.update(summary_rows=len(rates),aggregate_rows=len(agg),bootstrap_draws=400,summary_checks_passed=True)
    dump(V/'results/audit.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
