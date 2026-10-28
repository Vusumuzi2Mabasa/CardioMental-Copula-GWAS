"""Independent completeness, decision reconstruction and seed replay audit."""
from pathlib import Path
import json,csv,hashlib,subprocess
from datetime import datetime
import numpy as np
from benchmark import cells,codehash,evaluate,CFG,ROOT
from simulate import generate,EFFECTS
from model import METHODS


def read(path):
    with Path(path).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    frozen=subprocess.check_output(['git','-c','core.fsmonitor=false','show','fe4bb39:protocol/design.json'],cwd=ROOT)
    assert frozen==(ROOT/'protocol/design.json').read_bytes(),'Frozen design changed'
    for name in ['src/model.py','src/simulate.py','src/benchmark.py']:
        archived=subprocess.check_output(['git','-c','core.fsmonitor=false','show','fe4bb39:'+name],cwd=ROOT)
        assert archived==(ROOT/name).read_bytes(),'Frozen implementation changed: '+name
    identity=codehash();fingerprints=set();ncal=neval=nscore=neffects=0
    fails=cal_fails=panel_fails=0; details=[];replay=[]
    for c in cells():
        folder=ROOT/'results/cells'/f"cell_{c['id']:02d}"
        complete=json.loads((folder/'complete.json').read_text())
        assert complete['code_sha256']==identity
        for name,h in complete['outputs'].items():assert sha(folder/name)==h,(c['id'],name)
        freeze=json.loads((folder/'threshold_freeze.json').read_text())
        assert freeze['code_sha256']==identity
        assert freeze['calibration_sha256']==sha(folder/'calibration.csv')
        assert datetime.fromisoformat(freeze['frozen_before_evaluation_utc'])<datetime.fromisoformat(complete['completed_utc'])
        cal=read(folder/'calibration.csv');datasets=read(folder/'datasets.csv');rows=read(folder/'scores.csv')
        ff=read(folder/'null_fits.csv');ee=read(folder/'effects.csv')
        assert len(cal)==len(METHODS)*c['cal']
        expected=c['null']+4*c['alt']
        assert len(datasets)==expected and len(rows)==len(METHODS)*expected and len(ff)==3*expected
        assert len(ee)==5*c['effect_reps']*4
        for method in METHODS:
            cr=[r for r in cal if r['method']==method]
            assert {int(r['replicate']) for r in cr}==set(range(c['cal']))
            nfail=sum(r['status']!='ok' or not np.isfinite(float(r['p'])) for r in cr)
            assert nfail==freeze['calibration_failures'][method]
            cal_fails+=nfail
            cutoff=freeze['thresholds'][method]
            if nfail:assert cutoff is None
            else:assert cutoff==sorted(float(r['p']) for r in cr)[9]
            for scenario in EFFECTS:
                selected=[r for r in rows if r['method']==method and r['scenario']==scenario]
                count=c['null'] if scenario=='null' else c['alt']
                assert {int(r['replicate']) for r in selected}==set(range(count)) and len(selected)==count
                for r in selected:
                    p=float(r['p']);valid=r['status']=='ok'
                    assert np.isfinite(p)==valid
                    if valid:assert 0<=p<=1
                    assert int(r['nominal_05'])==int(p<.05)
                    assert int(r['nominal_01'])==int(p<.01)
                    if cutoff is not None:assert int(r['calibrated_05'])==int(p<cutoff)
                    else:assert r['calibrated_05']==''
                    ps=[float(r['p'+str(j)]) for j in range(CFG['markers'])]
                    assert int(r['panel_FWER_05'])==int(any(q<.05/CFG['markers'] for q in ps))
                    fails+=not valid;panel_fails+=r['panel_status']!='ok'
        for r in datasets:
            assert r['fingerprint'] not in fingerprints,'Duplicate evaluation dataset'
            fingerprints.add(r['fingerprint'])
        for r in ee:
            if r['status']=='ok':
                e=float(r['estimate']);t=float(r['truth'])
                for var,key in [('sandwich_se','sandwich_coverage'),('model_se','model_coverage')]:
                    se=float(r[var]);assert np.isfinite(se) and se>0
                    assert int(r[key])==int(abs(e-t)<=1.95996398454*se)
            assert int(r['replicate'])<c['effect_reps']
        # Regenerate a null and a shared-effect dataset in each stratum.
        for scenario in ['null','shared']:
            dat=generate(CFG['seed_evaluation'],c['id'],0,scenario,c,CFG['markers'])
            uid=hashlib.sha256(dat['g'].tobytes()+dat['y'].tobytes()+dat['d'].tobytes()).hexdigest()
            target=next(r for r in datasets if r['scenario']==scenario and int(r['replicate'])==0)
            assert uid==target['fingerprint']
            scores,_,_=evaluate(dat,False)
            differences=[]
            for method,z in scores.items():
                r=next(r for r in rows if r['scenario']==scenario and int(r['replicate'])==0 and r['method']==method)
                assert z['status']==r['status']
                if r['status']=='ok':
                    diff=max(abs(a-float(r['p'+str(j)])) for j,a in enumerate(z['p']))
                    assert diff<2e-6;differences.append(diff)
            replay.append(dict(cell=c['id'],scenario=scenario,max_p_error=max(differences,default=0)))
        ncal+=c['cal'];neval+=len(datasets);nscore+=len(rows);neffects+=len(ee)
        details.append(dict(cell=c['id'],evaluation_datasets=len(datasets),mean_case_fraction=np.mean([int(r['cases'])/c['n'] for r in datasets])))
    assert ncal==3200 and neval==9920 and nscore==69440 and neffects==6400
    result=dict(passed=True,protocol_commit='fe4bb394d4a33847d00d7559f2a1371b5342afad',code_sha256=identity,
        calibration_datasets=ncal,evaluation_datasets=neval,total_datasets=ncal+neval,
        evaluation_score_records=nscore,variant_scores=nscore*CFG['markers'],unique_evaluation_fingerprints=len(fingerprints),
        calibration_score_failures=cal_fails,evaluation_focal_score_failures=fails,evaluation_panel_score_failures=panel_fails,
        effect_interval_records=neffects,effect_fit_failures=0,
        replayed_datasets=len(replay),replay_max_p_error=max(r['max_p_error'] for r in replay),
        limitation='Independent R controls validate numerical calculations, not all statistical operating properties')
    result['effect_fit_failures']=sum(sum(r['status']!='ok' for r in read(p)) for p in (ROOT/'results/cells').glob('*/effects.csv'))//2
    (ROOT/'results/audit.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'results/replay_audit.json').write_text(json.dumps(replay,indent=2)+'\n')
    (ROOT/'results/dataset_audit.json').write_text(json.dumps(details,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
