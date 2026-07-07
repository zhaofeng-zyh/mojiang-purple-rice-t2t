source ~/miniconda3/etc/profile.d/conda.sh
WD="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/23_亚种归属_3KRG_A12"
REF="$WD/ref"; A12="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/A12"
NIP="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP/nip_genome.fa"
mkdir -p "$WD/proj"; cd "$WD/proj"
echo "== wait for bed $(date) =="
for i in $(seq 1 120); do
  sz=$(stat -c%s "$REF/pruned_v2.1.bed" 2>/dev/null || echo 0)
  [ "$sz" -ge 764000000 ] && { echo "bed ready $sz"; break; }
  sleep 15
done
# 1) ZN65 genotype VCF at pruned positions
conda activate cgsv
python3 - <<'PY'
import bisect
REF="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/23_亚种归属_3KRG_A12/ref"
A12="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/A12"
NIP="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP/nip_genome.fa"
pruned={}; bychr={}
for ln in open(REF+"/pruned_v2.1.bim"):
    p=ln.split(); c="Chr"+p[0]; pos=int(p[3]); pruned[(c,pos)]=(p[4].upper(),p[5].upper(),p[0]); bychr.setdefault(c,[]).append(pos)
for c in bychr: bychr[c].sort()
nipref={}; cur=None; seq=[]
def flush(c,s):
    if c in bychr:
        S="".join(s).upper()
        for pos in bychr[c]:
            if pos-1<len(S): nipref[(c,pos)]=S[pos-1]
for ln in open(NIP):
    if ln[0]==">": flush(cur,seq); cur=ln[1:].split()[0]; seq=[]
    else: seq.append(ln.strip())
flush(cur,seq)
# ZN65 var
Rint={}; V={}
for ln in open(A12+"/ZN65.var"):
    t=ln[:1]
    if t=="R":
        p=ln.split("\t"); Rint.setdefault(p[1],[]).append((int(p[2]),int(p[3])))
    elif t=="V":
        p=ln.rstrip("\n").split("\t")
        if len(p)>=8: 
            c=p[1]; pos=int(p[2])+1
            if (c,pos) in pruned and len(p[6])==1 and len(p[7])==1: V[(c,pos)]=p[7].upper()
for c in Rint: Rint[c].sort()
def callable_at(c,pos):
    iv=Rint.get(c); 
    if not iv: return False
    st=[a for a,b in iv]; i=bisect.bisect_right(st,pos-1)-1
    return i>=0 and iv[i][0]<pos<=iv[i][1]
o=open("zn65.vcf","w")
o.write("##fileformat=VCFv4.2\n")
for n in range(1,13): o.write(f"##contig=<ID={n}>\n")
o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tZN65\n")
nw=0
for (c,pos),(a1,a2,cn) in pruned.items():
    if (c,pos) not in nipref: continue
    if not callable_at(c,pos): continue
    ref=nipref[(c,pos)]
    if ref not in (a1,a2): continue
    alt=a2 if ref==a1 else a1
    zb=V.get((c,pos), ref)
    if zb==ref: gt="0/0"
    elif zb==alt: gt="1/1"
    else: continue
    o.write(f"{cn}\t{pos}\t{cn}:{pos}\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}\n")
    nw+=1
o.close(); print("ZN65 VCF records:",nw)
PY
conda activate plink
# 2) 3K -> clean bed with chr:pos IDs
plink2 --bed "$REF/pruned_v2.1.bed" --bim "$REF/pruned_v2.1.bim" --fam "$REF/pruned_v2.1.fam" \
  --set-all-var-ids @:# --make-bed --out k3 --allow-extra-chr >/dev/null 2>&1
# 3) ZN65 -> bed
plink2 --vcf zn65.vcf --set-all-var-ids @:# --make-bed --out zn65 --allow-extra-chr >/dev/null 2>&1
# 4) merge (intersect SNPs)
plink2 --bfile k3 --pmerge zn65 --make-bed --out merged --allow-extra-chr >/dev/null 2>&1 || \
plink --bfile k3 --bmerge zn65 --make-bed --out merged --allow-extra-chr >/dev/null 2>&1
echo "merged fam count:"; wc -l merged.fam 2>/dev/null
# 5) PCA
plink2 --bfile merged --pca 10 --out merged_pca --allow-extra-chr >/dev/null 2>&1
echo "=== ZN65 在 merged PCA 的坐标 ==="; grep -i zn65 merged_pca.eigenvec
echo "DONE $(date)"; touch "$WD/proj/PROJ_DONE"
