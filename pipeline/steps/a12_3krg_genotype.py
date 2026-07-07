import os,bisect,numpy as np
WD="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/23_亚种归属_3KRG_A12"
REF=WD+"/ref"; A12="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/A12"
NIP="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP/nip_genome.fa"
os.makedirs(WD+"/out",exist_ok=True)
# 1) pruned positions
pruned={}   # (chr,pos)->(A1,A2)
bychr={}
for ln in open(REF+"/pruned_v2.1.bim"):
    p=ln.split()
    c="Chr"+p[0]; pos=int(p[3]); a1=p[4].upper(); a2=p[5].upper()
    pruned[(c,pos)]=(a1,a2); bychr.setdefault(c,[]).append(pos)
for c in bychr: bychr[c].sort()
print("pruned SNPs:",len(pruned))
# 2) Nipponbare ref base at pruned positions
nipref={}
cur=None; seq=[]
def flush(c,s):
    if c and c in bychr:
        S="".join(s).upper()
        for pos in bychr[c]:
            if pos-1 < len(S): nipref[(c,pos)]=S[pos-1]
for ln in open(NIP):
    if ln[0]==">":
        flush(cur,seq); cur=ln[1:].split()[0]; seq=[]
    else: seq.append(ln.strip())
flush(cur,seq)
print("nip ref bases:",len(nipref))
samples=["ZN65","MH63","ZS97","CempoIreng","N22","Orufipogon"]
meta={"ZN65":"?/purple","MH63":"indica","ZS97":"indica","CempoIreng":"trop-japonica","N22":"aus","Orufipogon":"wild"}
# 3) genotype each sample at pruned positions
def load(samp):
    Rint={}; V={}
    for ln in open(f"{A12}/{samp}.var"):
        t=ln[:1]
        if t=="R":
            p=ln.split("\t"); Rint.setdefault(p[1],[]).append((int(p[2]),int(p[3])))
        elif t=="V":
            p=ln.rstrip("\n").split("\t")
            if len(p)<8: continue
            c=p[1]; pos=int(p[2])+1
            if (c,pos) in pruned:
                ref=p[6].upper(); alt=p[7].upper()
                if len(ref)==1 and len(alt)==1: V[(c,pos)]=alt
    for c in Rint: Rint[c].sort()
    return Rint,V
def callable_at(Rint,c,pos):
    iv=Rint.get(c); 
    if not iv: return False
    st=[a for a,b in iv]; i=bisect.bisect_right(st,pos-1)-1
    return i>=0 and iv[i][0]<pos<=iv[i][1]
geno={}   # sample-> dict pos->base
for s in samples:
    Rint,V=load(s); g={}
    for k in pruned:
        c,pos=k
        if k not in nipref: continue
        if callable_at(Rint,c,pos):
            g[k]=V.get(k, nipref[k])   # alt if variant else nip ref
    geno[s]=g
    print(f"  {s}: callable@pruned {len(g):,}")
# 4) keep positions callable in ALL samples + biallelic-consistent
order=["Nipponbare"]+samples
keep=[]
for k in pruned:
    if k not in nipref: continue
    if all(k in geno[s] for s in samples):
        a1,a2=pruned[k]; alleles={nipref[k]}|{geno[s][k] for s in samples}
        if alleles<= {a1,a2}:   # strictly biallelic at the 3K-defined alleles
            keep.append(k)
print("usable (callable in all 7, biallelic):",f"{len(keep):,}")
# 5) build alignment + 0/1 matrix vs Nipponbare
seqs={s:[] for s in order}; M=[]
for k in keep:
    rb=nipref[k]; seqs["Nipponbare"].append(rb); row=[0]
    for s in samples:
        b=geno[s][k]; seqs[s].append(b); row.append(0 if b==rb else 1)
    M.append(row)
open(WD+"/out/snp_3krg.fasta","w").write("".join(f">{s}\n{''.join(seqs[s])}\n" for s in order))
# pairwise p-dist to ZN65
print("\n=== ZN65 与各样本 p-距离 (3K curated markers) ===")
zi=order.index("ZN65")
res=[]
for s in order:
    if s=="ZN65": continue
    d=t=0
    for k in keep:
        a=seqs["ZN65"][keep.index(k)] if False else None
    # vectorized
G=np.array([[ (0 if seqs[s][i]==seqs["Nipponbare"][i] else 1) for i in range(len(keep))] for s in order])
znrow=G[zi]
for i,s in enumerate(order):
    if s=="ZN65": continue
    # p-distance on actual bases
    diff=sum(1 for j in range(len(keep)) if seqs["ZN65"][j]!=seqs[s][j])
    print(f"  ZN65 vs {s:12s} {meta.get(s,'jap-ref'):14s} {100.0*diff/len(keep):.3f}%  ({diff}/{len(keep)})")
# PCA
Mc=G.astype(float); Mc=Mc-Mc.mean(0)
U,S,Vt=np.linalg.svd(Mc,full_matrices=False); pc=U*S; ve=(S**2)/np.sum(S**2)*100
with open(WD+"/out/pca_3krg.tsv","w") as o:
    o.write("sample\tgroup\tPC1\tPC2\tPC3\n")
    for i,s in enumerate(order):
        o.write(f"{s}\t{meta.get(s,'jap-ref')}\t{pc[i,0]:.2f}\t{pc[i,1]:.2f}\t{pc[i,2]:.2f}\n")
print(f"\nPCA var: PC1 {ve[0]:.1f}% PC2 {ve[1]:.1f}% PC3 {ve[2]:.1f}%")
