import os,numpy as np,bisect
A12="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/A12"
os.chdir(A12)
samples=["ZN65","MH63","ZS97","CempoIreng","N22","Orufipogon"]
meta={"ZN65":"japonica/purple","MH63":"indica/white","ZS97":"indica/white","CempoIreng":"tropjap/black","N22":"aus/white","Orufipogon":"wild"}
Rreg={}; V={}; refbase={}; allpos=set()
for s in samples:
    reg={}; v={}
    for ln in open(s+".var"):
        t=ln[:1]
        if t=="R":
            p=ln.split("\t"); reg.setdefault(p[1],[]).append((int(p[2]),int(p[3])))
        elif t=="V":
            p=ln.rstrip("\n").split("\t")
            if len(p)<8: continue
            ref=p[6].upper(); alt=p[7].upper()
            if len(ref)==1 and len(alt)==1 and ref in "ACGT" and alt in "ACGT":
                k=(p[1],int(p[2])); v[k]=alt; allpos.add(k); refbase.setdefault(k,ref)
    for c in reg: reg[c].sort()
    Rreg[s]=reg; V[s]=v
# callable starts per chr for binary search
starts={s:{c:[a for a,b in reg] for c,reg in Rreg[s].items()} for s in samples}
def callable_in(s,c,pos):
    reg=Rreg[s].get(c); 
    if not reg: return False
    st=starts[s][c]; i=bisect.bisect_right(st,pos)-1
    return i>=0 and reg[i][0]<=pos<=reg[i][1]
# keep sites callable in ALL samples
keep=[]
for (c,pos) in allpos:
    if all(callable_in(s,c,pos) for s in samples): keep.append((c,pos))
print("union SNP:",f"{len(allpos):,}","  high-confidence (callable in ALL 6):",f"{len(keep):,}")
keep.sort()
MAX=200000
if len(keep)>MAX:
    step=max(1,len(keep)//MAX); keep=keep[::step][:MAX]
print("sites used:",f"{len(keep):,}")
order=["Nipponbare"]+samples
seqs={s:[] for s in order}; geno=[]
for k in keep:
    rb=refbase[k]; seqs["Nipponbare"].append(rb); row=[0]
    for s in samples:
        a=V[s].get(k,rb); seqs[s].append(a); row.append(0 if a==rb else 1)
    geno.append(row)
open("snp_aln_hc.fasta","w").write("".join(f">{s}\n{''.join(seqs[s])}\n" for s in order))
# per-sample HC alt counts
print("=== HC alt counts vs Nipponbare ===")
g=np.array(geno); 
for i,s in enumerate(samples): print(f"  {s:11s} {meta[s]:15s} {int(g[:,i+1].sum()):,}")
M=g.astype(float).T; M=M-M.mean(0)
U,S,Vt=np.linalg.svd(M,full_matrices=False); pc=U*S; ve=(S**2)/np.sum(S**2)*100
with open("pca_hc.tsv","w") as o:
    o.write("sample\tgroup\tPC1\tPC2\tPC3\n")
    for i,s in enumerate(order):
        o.write(f"{s}\t{('jap-ref' if s=='Nipponbare' else meta[s])}\t{pc[i,0]:.1f}\t{pc[i,1]:.1f}\t{pc[i,2]:.1f}\n")
print(f"PCA var: PC1 {ve[0]:.1f}% PC2 {ve[1]:.1f}% PC3 {ve[2]:.1f}%")
