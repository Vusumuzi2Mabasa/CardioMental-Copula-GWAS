"""Reproducible public Phase 3 haplotype extraction, without donor identifiers in outputs."""
from pathlib import Path
import argparse, csv, gzip, hashlib, json, urllib.request
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
V=ROOT/'v3'
BASE='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
VCF='ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'
PANEL='integrated_call_samples_v3.20130502.ALL.panel'
EXPECTED={VCF:'a90c16c4ff2b3196476d506ae13cb3047fae8670163c7c932c4b0239aef3daf5',PANEL:'b4023dc6ee2d62ee89c8d4d347db4d348e65518d66d346574cdae7a4bbd76858'}
WINDOWS=[(20000000,20100000),(30000000,30100000),(40000000,40100000)]
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def table(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--raw',type=Path,default=V/'work/raw');args=ap.parse_args()
    args.raw.mkdir(parents=True,exist_ok=True);sources=[]
    for file in [VCF,PANEL]:
        p=args.raw/file
        if not p.exists():urllib.request.urlretrieve(BASE+file,p)
        digest=sha(p)
        if digest!=EXPECTED[file]:raise ValueError('Public source checksum differs: '+file)
        sources.append(dict(file=file,url=BASE+file,sha256=digest,bytes=p.stat().st_size))
    pops=['YRI','LWK'];names=[];splits={};rng=np.random.default_rng(29103100)
    lines=[x.split() for x in (args.raw/PANEL).read_text().splitlines()[1:] if x.strip()]
    for pop in pops:
        ids=sorted(t[0] for t in lines if t[1]==pop);ids=list(rng.permutation(ids))
        splits[pop]={'reference':ids[:40],'study':ids[40:]};names.extend(ids)
    records=[[] for _ in WINDOWS];haps=[[] for _ in WINDOWS]
    with gzip.open(args.raw/VCF,'rt') as f:
        for line in f:
            if line.startswith('##'):continue
            if line.startswith('#CHROM'):
                allnames=line.rstrip().split('\t')[9:];ix=[9+allnames.index(n) for n in names];continue
            pos=int(line.split('\t',2)[1])
            if pos>WINDOWS[-1][1]:break
            wi=next((i for i,(a,b) in enumerate(WINDOWS) if a<=pos<=b),None)
            if wi is None:continue
            t=line.rstrip().split('\t');ref,alt=t[3:5]
            if t[6]!='PASS' or len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT' or {ref,alt} in ({'A','T'},{'C','G'}):continue
            gt=[t[k].split(':')[0] for k in ix]
            if any(x not in ('0|0','0|1','1|0','1|1') for x in gt):continue
            records[wi].append(dict(id=f'22:{pos}:{ref}:{alt}',position=pos,ref=ref,alt=alt))
            haps[wi].append([[int(g[0]),int(g[2])] for g in gt])
    (V/'reference').mkdir(parents=True,exist_ok=True)
    qc=[];vs=[];pairs=[]
    for wi in range(3):
        h=np.array(haps[wi],dtype=np.int8).transpose(1,2,0);blocks={}
        for pop in pops:
            for role,ids in splits[pop].items():blocks[pop+'_'+role]=h[[names.index(n) for n in ids]].reshape(-1,h.shape[2])
        good=np.ones(h.shape[2],bool)
        for b in blocks.values():good &= np.minimum(b.mean(0),1-b.mean(0))>=.10
        candidates=np.flatnonzero(good)
        if len(candidates)<64:raise ValueError('Fewer than 64 common markers')
        # Dense retained marker panel is independent of all simulated phenotypes.
        chosen=candidates[np.linspace(0,len(candidates)-1,64).round().astype(int)]
        arr={};metadata={}
        for pop in pops:
            for role in ['reference','study']:
                key=pop+'_'+role;b=blocks[key][:,chosen];arr[key]=b
                rr=np.corrcoef(b,rowvar=False);arr[key+'_R']=rr
                qc.append(dict(region=wi,population=pop,pool=role,donors=len(b)//2,markers=64,min_maf=float(np.minimum(b.mean(0),1-b.mean(0)).min()),max_r2=float(np.max((rr*rr-np.eye(64))))))
            # Causal pair selection uses external reference donors only, never outcomes.
            r2=arr[pop+'_reference_R']**2
            a,b=np.where(np.triu((r2>=.5)&(r2<.98),1))
            if len(a)==0:raise ValueError('No moderately/highly correlated causal pair')
            k=np.argmax(r2[a,b]);anchor,high=int(a[k]),int(b[k])
            low=int(np.argmin(np.where(np.arange(64)==anchor,np.inf,r2[anchor])))
            if r2[anchor,low]>.05:raise ValueError('No low LD control')
            metadata[pop]=dict(anchor=anchor,high=high,low=low,reference_r2_high=float(r2[anchor,high]),reference_r2_low=float(r2[anchor,low]),study_r2_high=float(arr[pop+'_study_R'][anchor,high]**2),study_r2_low=float(arr[pop+'_study_R'][anchor,low]**2))
            pairs.append(dict(region=wi,population=pop,**metadata[pop]))
        for j,i in enumerate(chosen):vs.append(dict(region=wi,index=j,**records[wi][i]))
        np.savez_compressed(V/f'reference/region_{wi}.npz',**arr)
        dump(V/f'reference/region_{wi}.json',metadata)
    table(V/'reference/qc.csv',qc);table(V/'reference/variants.csv',vs);table(V/'reference/causal_pairs.csv',pairs)
    dump(V/'reference/provenance.json',dict(build='GRCh37',collection='1000 Genomes Phase 3 v5b',sources=sources,windows=WINDOWS,split_seed=29103100,reference_donors_each=40,donor_overlap=0,donor_list_digests={p:{r:hashlib.sha256(('\n'.join(ns)+'\n').encode()).hexdigest() for r,ns in roles.items()} for p,roles in splits.items()},reuse='Public IGSR data; derived haplotypes have no sample identifiers. Cite IGSR and 1000 Genomes. No participant phenotypes used.'))
    print(json.dumps(pairs,indent=2))
if __name__=='__main__':main()
