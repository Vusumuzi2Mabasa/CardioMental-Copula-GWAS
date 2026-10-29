from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];sys.path.append(str(ROOT/'v3/src'))
from prepare import sha,dump
V=ROOT/'v4'
def inputs():
    paths=['v4/protocol/design.json','v3/protocol/source_lock.json']
    paths += ['v4/src/'+p for p in ['kernel.py','reml.py','polygenic.py','export_finemap.py','finemap.R','finemap_summary.py']]
    paths += [str(p.relative_to(ROOT)).replace('\\','/') for p in sorted((V/'reference').glob('*')) if p.is_file()]
    return paths
def identity():return hashlib.sha256(json.dumps({p:sha(ROOT/p) for p in inputs()},sort_keys=True).encode()).hexdigest()
def verify():
    lock=json.loads((V/'protocol/source_lock.json').read_text());assert identity()==lock['identity'],'Frozen v4 inputs changed';return lock
def main():
    path=V/'protocol/source_lock.json'
    if path.exists():verify();print('Version four freeze verified');return
    import subprocess
    dump(path,dict(identity=identity(),sha256={p:sha(ROOT/p) for p in inputs()},md5={p:hashlib.md5((ROOT/p).read_bytes()).hexdigest() for p in inputs()},parent_commit=subprocess.check_output(['git','-c','core.fsmonitor=false','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),stage='Before version four fine-mapping and polygenic experiments'))
    print(identity())
if __name__=='__main__':main()
