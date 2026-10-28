"""Freeze tested source and protocol before the new evaluation starts."""
import json
from run import ROOT,REPO,CORE,sha,identity

def main():
    files=CORE+['v2/src/test_oracle.py','v2/src/crosscheck.R','v2/src/pilot.py','v2/results/pilot.csv','v2/results/pilot_selection.json']
    for name in ['numerical_controls.json','R_validation.json']:
        assert json.loads((ROOT/'results'/name).read_text())['passed']
    lock=dict(code_sha256=identity(),files={p:sha(REPO/p) for p in files})
    out=ROOT/'protocol/source_lock.json'
    if out.exists():assert json.loads(out.read_text())==lock,'Do not replace a frozen protocol'
    else:out.write_text(json.dumps(lock,indent=2)+'\n')
    print('Frozen',lock['code_sha256'])
if __name__=='__main__':main()
