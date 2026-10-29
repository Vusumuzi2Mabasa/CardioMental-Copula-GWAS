"""Freeze the design before any main calibration or evaluation data are generated."""
from prepare import ROOT,V,sha,dump
from run import identity
import subprocess
def main():
    dest=V/'protocol/source_lock.json'
    if dest.exists():
        import json
        if json.loads(dest.read_text())['identity']!=identity():raise ValueError('Frozen source changed')
        print('Existing freeze verified');return
    paths=[str(p.relative_to(ROOT)).replace('\\','/') for p in (V/'src').glob('*.py')]
    dump(dest,dict(identity=identity(),sources={p:sha(ROOT/p) for p in paths},parent_commit=subprocess.check_output(['git','-c','core.fsmonitor=false','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),stage='Before main simulations'))
    print(identity())
if __name__=='__main__':main()
