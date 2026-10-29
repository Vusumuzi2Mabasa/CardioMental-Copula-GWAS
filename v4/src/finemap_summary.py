"""Freeze level-two calibration before its evaluation, then audit all retained fits."""
from pathlib import Path
import csv,json,argparse,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
import sys
sys.path.insert(0,str(ROOT/'v3/src'))
from prepare import dump,table,sha
from model import SCENARIOS,NO_SHARED
V=ROOT/'v4'
def load(phase):
    rows=[]
    for p in sorted((V/'results/finemap_chunks').glob(phase+'_*.csv')):
        if p.stem.endswith('_sets'):continue
        with p.open() as f:
            for r in csv.DictReader(f):
                for k in ['region','replicate','max_pp4','signal_pairs','c_sets','b_sets','c_covered','b_covered','c_false_sets','b_false_sets','c_max_cs','b_max_cs']:r[k]=float(r[k])
                rows.append(r)
    return rows
def tag(r):return f"{r['population']}_r{int(r['region'])}_{r['regime']}_{r['LD']}"
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--freeze',action='store_true');args=ap.parse_args()
    cal=load('cal');assert len(cal)==2400
    ts={}
    for key in sorted(set(map(tag,cal))):
        a=[r for r in cal if tag(r)==key];assert len(a)==100
        failures=sum(r['status']!='ok' for r in a)
        ts[key]=dict(active=failures==0,failures=failures,cutoff=float(max(r['max_pp4'] for r in a)),calibration_per_architecture=20)
    path=V/'results/finemap_thresholds.json'
    if args.freeze:
        value=dict(thresholds=ts,alpha=.05,rank=20,strict=True,averaged_bound=1/21,calibration_rows=2400)
        if path.exists():assert json.loads(path.read_text())==value
        else:dump(path,value)
        print('Fine-mapping thresholds frozen',len(ts),'strata');return
    assert json.loads(path.read_text())['thresholds']==ts
    rows=load('eval');assert len(rows)==3840
    summaries=[]
    for key in ts:
        for s in SCENARIOS:
            a=[r for r in rows if tag(r)==key and r['scenario']==s];assert len(a)==20
            t=ts[key];first=a[0]
            ctrue=0 if s in ['null','mental_only'] else 1
            btrue=0 if s in ['null','cardio_only'] else (2 if s=='shared_plus_distinct' else 1)
            summaries.append(dict(population=first['population'],region=int(first['region']),regime=first['regime'],LD=first['LD'],scenario=s,attempted=20,failures=sum(r['status']!='ok' for r in a),calibration_active=t['active'],conventional_rate=sum(r['status']=='ok' and r['max_pp4']>.8 for r in a)/20,calibrated_rate=sum(t['active'] and r['status']=='ok' and r['max_pp4']>t['cutoff'] for r in a)/20,c_credible_coverage=sum(r['c_covered'] for r in a)/(20*ctrue) if ctrue else 0,b_credible_coverage=sum(r['b_covered'] for r in a)/(20*btrue) if btrue else 0,c_false_sets=sum(r['c_false_sets'] for r in a),b_false_sets=sum(r['b_false_sets'] for r in a),c_mean_max_set_size=float(np.mean([r['c_max_cs'] for r in a])),b_mean_max_set_size=float(np.mean([r['b_max_cs'] for r in a]))))
    table(V/'results/finemap_rates.csv',summaries)
    aggregate=[]
    for regime in ['ideal','realistic']:
        for ld in ['matched_reference','in_sample']:
            for s in SCENARIOS:
                a=[r for r in summaries if r['regime']==regime and r['LD']==ld and r['scenario']==s]
                aggregate.append(dict(regime=regime,LD=ld,scenario=s,attempted=120,failures=sum(r['failures'] for r in a),inactive_calibrations=sum(not r['calibration_active'] for r in a),**{k:float(np.mean([r[k] for r in a])) for k in ['conventional_rate','calibrated_rate','c_credible_coverage','b_credible_coverage','c_mean_max_set_size','b_mean_max_set_size']}))
    table(V/'results/finemap_aggregate.csv',aggregate)
    dump(V/'results/finemap_audit.json',dict(passed=True,calibration_analyses=len(cal),evaluation_analyses=len(rows),distinct_calibration_datasets=1200,distinct_evaluation_datasets=1920,LD_arms=2,calibration_failures=sum(r['status']!='ok' for r in cal),evaluation_failures=sum(r['status']!='ok' for r in rows),inactive_calibrations=sum(not t['active'] for t in ts.values()),not_independent_of_v3='Fixed first 20 replicates of the original version three calibration/evaluation datasets; no pooling as new independent replicates',threshold_sha256=sha(path)))
    print('Fine-mapping summaries audited')
if __name__=='__main__':main()
