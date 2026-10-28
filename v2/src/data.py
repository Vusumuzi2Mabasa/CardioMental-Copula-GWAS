"""Memory-bounded loading with strict replicate and method identity checks."""
from pathlib import Path
import numpy as np
from run import ROOT,config,cells,scenarios,read
from oracle import METHODS

def load(cell,phase,scenario):
    cfg=config();sid=0 if phase=='calibration' else [s['name'] for s in scenarios(cell)].index(scenario)
    n=cfg['calibration_null'] if phase=='calibration' else next(s['nrep'] for s in scenarios(cell) if s['name']==scenario)
    p=np.full((n,len(METHODS)),np.nan);seconds=np.zeros_like(p);seen=np.zeros_like(p,dtype=bool)
    good=np.zeros_like(p,dtype=bool);fingerprints=['']*n;cases=np.zeros(n,dtype=int)
    for path in sorted((ROOT/'results/chunks'/f"cell_{cell['id']}").glob(f'{phase}_{sid}_*.csv')):
        for r in read(path):
            i=int(r['replicate']);j=METHODS.index(r['method'])
            assert 0<=i<n and not seen[i,j],('duplicate/out of range',path,i,j)
            assert int(r['cell'])==cell['id'] and r['phase']==phase and r['scenario']==scenario
            seen[i,j]=True;p[i,j]=float(r['p']);seconds[i,j]=float(r['seconds'])
            good[i,j]=r['status']=='ok' and np.isfinite(p[i,j]) and 0<=p[i,j]<=1
            assert not fingerprints[i] or fingerprints[i]==r['fingerprint']
            fingerprints[i]=r['fingerprint'];cases[i]=int(r['cases'])
    assert seen.all(),('Missing records',cell,phase,scenario,int(seen.sum()))
    return dict(p=p,seconds=seconds,good=good,fingerprints=fingerprints,cases=cases,n=n)
