"""Frozen-source and full output audit with independent cutoffs and seeded replay."""
import json,zipfile,hashlib,subprocess
from datetime import datetime
import numpy as np
from run import ROOT,REPO,CORE,config,cells,scenarios,sha,identity,read
from data import load
from oracle import METHODS,evaluate
from generate import generate

def run():
    cfg=config();lock=json.loads((ROOT/'protocol/source_lock.json').read_text())
    assert identity()==lock['code_sha256']
    commit=subprocess.check_output(['git','-c','core.fsmonitor=false','log','--format=%H','-1','--','v2/protocol/source_lock.json'],cwd=REPO,text=True).strip()
    for path,h in lock['files'].items():
        assert sha(REPO/path)==h,path
        saved=subprocess.check_output(['git','-c','core.fsmonitor=false','show',commit+':'+path],cwd=REPO)
        assert hashlib.sha256(saved).hexdigest()==h,path
    # Verify every byte of the preserved complete version-one archive.
    with zipfile.ZipFile(REPO/'baseline/CardioMental_v1_48aa582.zip') as z:
        original=z.read('checksums_public.sha256').decode()
        for line in original.splitlines():
            expected,name=line.split('  ',1)
            assert hashlib.sha256(z.read(name)).hexdigest()==expected,name
            if name.startswith(('src/','results/','protocol/','figures/','report/')):
                assert sha(REPO/name)==expected,('Original file changed',name)
        assert z.read('README.md')==(REPO/'baseline/README_v1.md').read_bytes()
        assert z.read('report/CardioMental_Copula_GWAS_Report.docx')==(REPO/'baseline/report/CardioMental_Copula_GWAS_Report_v1.docx').read_bytes()
    freeze=json.loads((ROOT/'results/thresholds.json').read_text())
    assert freeze['code_sha256']==identity()
    timefreeze=datetime.fromisoformat(freeze['frozen_before_evaluation_utc'])
    for p,h in freeze['calibration_hashes'].items():assert sha(ROOT/p)==h,p
    counts=dict(calibration_datasets=0,evaluation_datasets=0,method_records=0,failures=0)
    allids=set();replays=[];reports=[]
    for cell in cells():
        for done in sorted((ROOT/'results/chunks'/f"cell_{cell['id']}").glob('*.json')):
            d=json.loads(done.read_text());csv=done.with_suffix('.csv')
            assert d['code_sha256']==identity() and d['csv_sha256']==sha(csv)
            rr=read(csv);assert len(rr)==d['rows']==(d['stop']-d['start'])*len(METHODS)
            assert sum(r['status']!='ok' for r in rr)==d['failures']
            if d['phase']=='evaluation':assert datetime.fromisoformat(d['completed_utc'])>timefreeze
        cal=load(cell,'calibration','null');counts['calibration_datasets']+=cal['n']
        for j,m in enumerate(METHODS):
            for alpha in cfg['alpha_grid']:
                r=next(t for t in freeze['thresholds'] if t['cell']==cell['id'] and t['method']==m and t['alpha']==alpha)
                rank=max(1,int(np.floor(alpha*(cal['n']+1))))
                fails=int((~cal['good'][:,j]).sum())
                assert r['rank']==rank and r['failures']==fails
                assert r['cutoff']==(None if fails else float(np.sort(cal['p'][:,j])[rank-1]))
        for phase,ss in [('calibration',[dict(name='null',beta=[0.,0.])]),('evaluation',scenarios(cell))]:
            for sid,sc in enumerate(ss):
                dat=cal if phase=='calibration' else load(cell,phase,sc['name'])
                for uid in dat['fingerprints']:assert uid not in allids;allids.add(uid)
                counts['method_records']+=dat['n']*len(METHODS);counts['failures']+=int((~dat['good']).sum())
                if phase=='evaluation':counts['evaluation_datasets']+=dat['n']
                reports.append(dict(**cell,phase=phase,scenario=sc['name'],datasets=dat['n'],failures=int((~dat['good']).sum()),
                                    mean_case_fraction=float(dat['cases'].mean()/cell['n'])))
                # Replay first replicate of every stream; all statuses compared.
                sample=generate(cfg['seed_calibration'] if phase=='calibration' else cfg['seed_evaluation'],cell,0,sc['beta'],stream=sid)
                from run import fingerprint
                assert fingerprint(sample)==dat['fingerprints'][0]
                again=evaluate(sample);diff=[]
                for j,m in enumerate(METHODS):
                    assert (again[m]['status']=='ok')==bool(dat['good'][0,j])
                    if dat['good'][0,j]:diff.append(abs(again[m]['p']-dat['p'][0,j]))
                assert max(diff,default=0)<1e-10
                replays.append(dict(cell=cell['id'],phase=phase,scenario=sc['name'],max_p_error=max(diff,default=0)))
    assert counts['calibration_datasets']==20000 and counts['evaluation_datasets']==72000
    assert counts['method_records']==828000 and len(allids)==92000
    result=dict(passed=True,protocol_commit=commit,code_sha256=identity(),**counts,total_datasets=92000,
                pilot_datasets=640,unique_dataset_fingerprints=len(allids),replayed_datasets=len(replays),
                replay_max_p_error=max(r['max_p_error'] for r in replays),baseline_bytes_preserved=True,
                calibration_frozen_before_evaluation=True,numerical_controls=json.loads((ROOT/'results/numerical_controls.json').read_text()),
                independent_R=json.loads((ROOT/'results/R_validation.json').read_text()))
    (ROOT/'results/audit.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'results/replay_audit.json').write_text(json.dumps(replays,indent=2)+'\n')
    (ROOT/'results/dataset_summary.json').write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':run()
