import glob,os,numpy as np
A12="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/A12"
os.chdir(A12)
samples=["ZN65","MH63","ZS97","CempoIreng","N22","Orufipogon"]  # Nipponbare=reference(all ref)
meta={"ZN65":"japonica/purple","MH63":"indica/white","ZS97":"indica/white",
      "CempoIreng":"tropjaponica/black","N22":"aus/white","Orufipogon":"wild"}
data={}; refbase={}; allpos=set(); counts={}
for s in samples:
    d={}; n=0
    for ln in open(s+".var"):
        if ln[:1]!="V": continue
        p=ln.rstrip("\n").split("\t")
        if len(p)<8: continue
        ref=p[6].upper(); alt=p[7].upper()
        if len(ref)==1 and len(alt)==1 and ref in "ACGT" and alt in "ACGT":
            key=(p[1],p[2]); d[key]=alt; allpos.add(key); n+=1
            if key not in refbase: refbase[key]=ref
    data[s]=d; counts[s]=n
print("=== per-sample SNPs vs Nipponbare (raw divergence) ===")
for s in samples: print(f"  {s:12s} {meta[s]:18s} {counts[s]:,}")
print("union SNP sites:", f"{len(allpos):,}")

# keep parsimony-useful sites; subsample deterministically for tractability
pos=sorted(allpos)
MAX=300000
if len(pos)>MAX:
    step=max(1,len(pos)//MAX); pos=pos[::step][:MAX]
print("sites used for tree/PCA:", f"{len(pos):,}")

order=["Nipponbare"]+samples
seqs={s:[] for s in order}
geno=[]  # rows=sites, cols=taxa (0 ref / 1 alt)
for k in pos:
    rb=refbase[k]; seqs["Nipponbare"].append(rb); row=[0]
    for s in samples:
        a=data[s].get(k,rb); seqs[s].append(a); row.append(0 if a==rb else 1)
    geno.append(row)
with open("snp_aln.fasta","w") as o:
    for s in order: o.write(">"+s+"\n"+"".join(seqs[s])+"\n")
print("snp_aln.fasta written:", sum(len(seqs[s]) for s in order), "chars")

# PCA (taxa x sites)
M=np.array(geno,dtype=float).T
M=M-M.mean(0)
U,S,Vt=np.linalg.svd(M,full_matrices=False)
pc=U*S
varexp=(S**2)/np.sum(S**2)*100
with open("pca.tsv","w") as o:
    o.write("sample\tgroup\tPC1\tPC2\tPC3\n")
    for i,s in enumerate(order):
        g=("japonica-ref" if s=="Nipponbare" else meta[s])
        o.write(f"{s}\t{g}\t{pc[i,0]:.2f}\t{pc[i,1]:.2f}\t{pc[i,2]:.2f}\n")
print(f"PCA var explained: PC1 {varexp[0]:.1f}% PC2 {varexp[1]:.1f}% PC3 {varexp[2]:.1f}%")
