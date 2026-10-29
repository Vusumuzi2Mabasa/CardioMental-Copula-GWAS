"""Public LD templates and independent synthetic polygenic relationship kernels."""
from pathlib import Path
import sys,argparse,json,gzip,hashlib
import numpy as np
from scipy.linalg import eigh
ROOT=Path(__file__).resolve().parents[2];sys.path.append(str(ROOT/'v3/src'))
from prepare import VCF,PANEL,EXPECTED,sha,dump,table
V=ROOT/'v4'
def reference(raw):
    for fn in [VCF,PANEL]:assert sha(raw/fn)==EXPECTED[fn]
    panel=[x.split() for x in (raw/PANEL).read_text().splitlines()[1:] if x.strip()]
    rng=np.random.default_rng(29103100);pops=['YRI','LWK'];names=[];by_pop={}
    for pop in pops:
        ids=sorted(t[0] for t in panel if t[1]==pop);ids=list(rng.permutation(ids))[40:]
        by_pop[pop]=ids;names.extend(ids)
    windows=[(20000000+1000000*k,20100000+1000000*k) for k in range(22)]
    hs=[[] for _ in windows];positions=[[] for _ in windows]
    with gzip.open(raw/VCF,'rt') as f:
        for line in f:
            if line.startswith('##'):continue
            if line.startswith('#CHROM'):
                s=line.rstrip().split('\t')[9:];ix=[9+s.index(n) for n in names];continue
            pos=int(line.split('\t',2)[1])
            if pos>windows[-1][1]:break
            k=(pos-20000000)//1000000
            if not 0<=k<22 or not windows[k][0]<=pos<=windows[k][1]:continue
            t=line.rstrip().split('\t');ref,alt=t[3:5]
            if t[6]!='PASS' or len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT' or {ref,alt} in ({'A','T'},{'C','G'}):continue
            gt=[t[j].split(':')[0] for j in ix]
            if any(g not in ('0|0','0|1','1|0','1|1') for g in gt):continue
            hs[k].append([[int(g[0]),int(g[2])] for g in gt]);positions[k].append(dict(position=pos,ref=ref,alt=alt))
    arrays={};variants=[];qc=[]
    for k in range(22):
        h=np.array(hs[k],dtype=np.int8).transpose(1,2,0);blocks={};good=np.ones(h.shape[2],bool)
        for pop in pops:
            b=h[[names.index(x) for x in by_pop[pop]]].reshape(-1,h.shape[2]);blocks[pop]=b
            good &= np.minimum(b.mean(0),1-b.mean(0))>=.10
        available=np.flatnonzero(good);assert len(available)>=32
        idx=available[np.linspace(0,len(available)-1,32).round().astype(int)]
        for pop,b in blocks.items():
            arrays[f'{pop}_{k}']=b[:,idx]
            qc.append(dict(population=pop,block=k,donors=len(b)//2,markers=32,min_maf=float(np.minimum(b[:,idx].mean(0),1-b[:,idx].mean(0)).min())))
        for j,i in enumerate(idx):variants.append(dict(synthetic_block=k,index=j,source_chromosome=22,**positions[k][i]))
    (V/'reference').mkdir(parents=True,exist_ok=True)
    np.savez_compressed(V/'reference/polygenic_haplotypes.npz',**arrays)
    table(V/'reference/polygenic_variants.csv',variants);table(V/'reference/polygenic_qc.csv',qc)
    dump(V/'reference/polygenic_provenance.json',dict(source='Same checksum-verified public Phase 3 GRCh37 VCF and disjoint study donor pools as version three',windows=windows,synthetic_blocks=22,markers_per_block=32,source_chromosome=22,important='Independent sampling across these windows defines a synthetic polygenic genome. These are not real chromosome 1 through 22 data; no real whole-genome African LD validation is claimed.'))
def kernel(pop,n):
    path=V/f'work/kernels/{pop}_n{n}.npz'
    if path.exists():return path
    rng=np.random.default_rng(np.random.SeedSequence([29103110,['YRI','LWK'].index(pop),n]));panels=np.load(V/'reference/polygenic_haplotypes.npz');blocks=[]
    for k in range(22):
        h=panels[f'{pop}_{k}'];g=h[rng.integers(len(h),size=n)]+h[rng.integers(len(h),size=n)]
        f=h.mean(0);blocks.append((g-2*f)/np.sqrt(2*f*(1-f)))
    Z=np.column_stack(blocks);Z=Z/np.sqrt(np.mean(Z*Z,axis=1))[:,None]
    K=Z@Z.T/Z.shape[1];assert max(abs(np.diag(K)-1))<1e-12
    X=np.column_stack([np.ones(n),rng.normal(size=n)])
    Q=np.linalg.qr(X,mode='complete')[0][:,2:];lam,U=eigh(Q.T@K@Q,check_finite=True)
    assert lam.min()>-1e-8
    lam=np.maximum(lam,0) # Numerical negative zero only, after explicit PSD check.
    B=Q@U;path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,Z=Z,B=B,eigenvalues=lam)
    public=V/f'reference/kernel_{pop}_n{n}.csv'
    if public.exists():
        import csv
        with public.open() as f:old=np.array([float(r['eigenvalue']) for r in csv.DictReader(f)])
        assert np.allclose(old,lam,atol=1e-9,rtol=1e-10),'Rebuilt kernel spectrum differs'
    else:table(public,[dict(mode=i,eigenvalue=float(v)) for i,v in enumerate(lam)])
    print('kernel prepared',pop,n,flush=True);return path
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--raw',type=Path);args=ap.parse_args()
    if args.raw:reference(args.raw)
    for pop in ['YRI','LWK']:
        for n in [600,1200]:kernel(pop,n)
if __name__=='__main__':main()
