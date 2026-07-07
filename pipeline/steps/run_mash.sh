source ~/miniconda3/etc/profile.d/conda.sh; conda activate mash
A12="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/A12"
GEN="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/genomes"
B="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定"
cd "$A12"
ZN65="/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa"
NIP="$B/08_全基因组SV_ZN65vsNIP/nip_genome.fa"
MH63="$B/10_多基因组Kala4/MH63RS3.fasta"; ZS97="$B/10_多基因组Kala4/ZS97RS3.fasta"
# label copies via symlink names so mash uses clean ids
ln -sf "$ZN65" ZN65.fa; ln -sf "$NIP" Nipponbare.fa; ln -sf "$MH63" MH63.fa; ln -sf "$ZS97" ZS97.fa
ln -sf "$GEN/CempoIreng_black.fna" CempoIreng.fa; ln -sf "$GEN/N22_aus.fna" N22.fa; ln -sf "$GEN/Orufipogon_wild.fna" Orufipogon.fa
mash sketch -k 21 -s 100000 -o panel ZN65.fa Nipponbare.fa MH63.fa ZS97.fa CempoIreng.fa N22.fa Orufipogon.fa 2>/dev/null
mash dist panel.msh panel.msh > mashdist.tsv 2>/dev/null
echo "=== mash 距离矩阵 (越小越近) ==="
python3 - <<'PY'
import numpy as np
names=[]; idx={}
rows={}
for ln in open("mashdist.tsv"):
    a,b,d,p,sh=ln.split("\t")
    a=a.replace(".fa",""); b=b.replace(".fa","")
    rows[(a,b)]=float(d)
    if a not in idx: idx[a]=len(names); names.append(a)
n=len(names); D=np.zeros((n,n))
for (a,b),d in rows.items(): D[idx[a],idx[b]]=d
# print matrix
print("\t"+"\t".join(f"{x[:6]}" for x in names))
for i,a in enumerate(names):
    print(f"{a[:11]:11s}\t"+"\t".join(f"{D[i,j]*100:.3f}" for j in range(n)))
# distance to ZN65 and to Nipponbare, sorted
zi=idx["ZN65"]; ni=idx["Nipponbare"]
print("\n=== 各基因组到 ZN65 的 mash 距离(%) 升序 ===")
for j in sorted(range(n),key=lambda j:D[zi,j]):
    if j!=zi: print(f"  {names[j]:12s} {D[zi,j]*100:.3f}")
# NJ tree (manual)
import copy
def nj(D,names):
    D=D.astype(float).copy(); nodes=list(names); 
    while len(nodes)>2:
        m=len(nodes); r=D.sum(1)
        Q=np.full((m,m),np.inf)
        for i in range(m):
            for j in range(m):
                if i!=j: Q[i,j]=(m-2)*D[i,j]-r[i]-r[j]
        i,j=np.unravel_index(np.argmin(Q),Q.shape)
        if i>j: i,j=j,i
        di=0.5*D[i,j]+ (r[i]-r[j])/(2*(m-2)) if m>2 else 0.5*D[i,j]
        dj=D[i,j]-di
        new=f"({nodes[i]}:{max(di,0):.4f},{nodes[j]}:{max(dj,0):.4f})"
        nd=np.array([0.5*(D[i,k]+D[j,k]-D[i,j]) for k in range(m)])
        keep=[k for k in range(m) if k!=i and k!=j]
        Dn=np.zeros((len(keep)+1,len(keep)+1)); nn=[nodes[k] for k in keep]+[new]
        for a,ka in enumerate(keep):
            for b,kb in enumerate(keep): Dn[a,b]=D[ka,kb]
            Dn[a,-1]=Dn[-1,a]=nd[ka]
        D=Dn; nodes=nn
    return f"({nodes[0]}:{D[0,1]/2:.4f},{nodes[1]}:{D[0,1]/2:.4f});"
print("\n=== NJ 树 (mash 距离) ===")
print(nj(D,names))
PY
