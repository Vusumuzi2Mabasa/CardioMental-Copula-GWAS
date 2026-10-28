"""Independent pilot chooses two scales by a frozen method-neutral rule."""
from pathlib import Path
import json,csv,time
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from oracle import evaluate
from generate import generate

ROOT=Path(__file__).resolve().parents[1]
ARCH={'continuous_only':(.2,0.),'binary_only':(0.,.35),'shared':(.2,.35),'opposing':(.2,-.35)}
def block(env,scale):
    cell=dict(id=0 if env=='constant' else 1,n=800,maf=.25,prevalence=.15,environment=env)
    rows=[]
    for j,(arch,b) in enumerate(ARCH.items()):
        for rep in range(20):
            dat=generate(29102001,cell,rep,tuple(scale*z for z in b),stream=j)
            ss=evaluate(dat)
            for m,z in ss.items():rows.append(dict(environment=env,scale=scale,architecture=arch,replicate=rep,method=m,**z))
    return rows

def main():
    rows=[];began=time.perf_counter()
    with ProcessPoolExecutor(max_workers=2) as pool:
        fut=[pool.submit(block,e,s) for e in ['constant','varying'] for s in [.25,.4,.55,.7]]
        for f in as_completed(fut):rows+=f.result();print('pilot rows',len(rows),flush=True)
    with (ROOT/'results/pilot.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    # Only the simpler prespecified marginal benchmark selects the signal range.
    rates={s:np.mean([r['p']<.05 for r in rows if r['scale']==s and r['method']=='Marginal covariance-aware']) for s in [.25,.4,.55,.7]}
    candidates=[(a,b) for a in rates for b in rates if a<b]
    selected=min(candidates,key=lambda ab:((rates[ab[0]]-.30)**2+(rates[ab[1]]-.65)**2,ab))
    result=dict(seed=29102001,replicates_per_cell=20,pooled_marginal_power=rates,selected_scales=selected,
                rule='Distinct ascending scales minimizing squared distance of pooled marginal nominal power to 0.30 and 0.65; ties lexicographic',
                datasets=len(rows)//9,seconds=time.perf_counter()-began)
    (ROOT/'results/pilot_selection.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
if __name__=='__main__':main()
