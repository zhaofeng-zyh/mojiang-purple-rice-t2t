source ~/miniconda3/etc/profile.d/conda.sh
WD="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/23_亚种归属_3KRG_A12"
REF="$WD/ref"; P="$WD/proj"; cd "$P"
echo "== wait clean bed $(date) =="
for i in $(seq 1 160); do
  sz=$(stat -c%s "$REF/pruned_v2.1.bed" 2>/dev/null || echo 0)
  [ "$sz" = "764770359" ] && { echo "bed clean $sz"; break; }
  sleep 15
done
conda activate plink
plink2 --bed "$REF/pruned_v2.1.bed" --bim "$REF/pruned_v2.1.bim" --fam "$REF/pruned_v2.1.fam" \
  --set-all-var-ids @:# --rm-dup force-first --make-bed --out k3 --allow-extra-chr >/dev/null 2>&1
echo "k3:"; wc -l k3.fam k3.bim 2>/dev/null
plink2 --bfile zn65 --rm-dup force-first --make-bed --out zn65b --allow-extra-chr >/dev/null 2>&1
# common variants then bmerge (plink1.9)
plink --bfile k3 --bmerge zn65b --make-bed --out merged --allow-extra-chr >/dev/null 2>merr.txt
if [ ! -f merged.fam ]; then
  echo "bmerge issue, trying with missnp exclude:"; tail -3 merr.txt
  if [ -f merged-merge.missnp ]; then
    plink --bfile k3 --exclude merged-merge.missnp --make-bed --out k3f --allow-extra-chr >/dev/null 2>&1
    plink --bfile zn65b --exclude merged-merge.missnp --make-bed --out zn65f --allow-extra-chr >/dev/null 2>&1
    plink --bfile k3f --bmerge zn65f --make-bed --out merged --allow-extra-chr >/dev/null 2>&1
  fi
fi
echo "merged:"; wc -l merged.fam 2>/dev/null
plink2 --bfile merged --pca 10 --out merged_pca --allow-extra-chr >/dev/null 2>&1
conda activate cgsv
python3 - <<'PY'
import numpy as np
# eigenvec
ev={}
for ln in open("merged_pca.eigenvec"):
    p=ln.split()
    if p[0]=="FID" or p[0]=="#FID": continue
    ev[p[1]]=np.array([float(x) for x in p[2:7]])
# Qmatrix subpop (argmax)
import csv
hdr=None; sub={}
for row in csv.reader(open("/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/23_亚种归属_3KRG_A12/ref/Qmatrix-k9-3kRG.csv")):
    if hdr is None: hdr=row[1:]; continue
    acc=row[0]; vals=[float(x) for x in row[1:]]
    sub[acc]=hdr[int(np.argmax(vals))]
groups={"ind1A":"indica","ind1B":"indica","ind2":"indica","ind3":"indica","temp":"japonica(temperate)","trop1":"japonica(tropical)","trop2":"japonica(tropical)","aus":"aus","aro":"aromatic"}
if "ZN65" not in ev:
    print("ZN65 not in PCA!"); print(list(ev)[:3]); raise SystemExit
z=ev["ZN65"]
# nearest 50 accessions
d=[]
for a,v in ev.items():
    if a=="ZN65": continue
    if a in sub: d.append((np.sqrt(((v-z)**2).sum()),a,sub[a]))
d.sort()
from collections import Counter
near=[g for _,_,g in d[:30]]
c=Counter([groups.get(g,g) for g in near])
print("ZN65 PC1-5:", np.round(z,3))
print("\n=== ZN65 最近 30 个 3K accession 的亚种 ===")
for di,a,g in d[:12]: print(f"  {a}  {groups.get(g,g):22s} dist={di:.4f}")
print("\n=== 最近 30 的亚种计票 ===")
for k,v in c.most_common(): print(f"  {k}: {v}")
# subpop centroids
cen={}
for a,v in ev.items():
    if a in sub:
        g=groups.get(sub[a],sub[a]); cen.setdefault(g,[]).append(v)
print("\n=== ZN65 到各亚种质心距离 (升序) ===")
cd=sorted(((np.sqrt(((np.mean(np.array(vs),0)-z)**2).sum()),g) for g,vs in cen.items()))
for di,g in cd: print(f"  {g:22s} {di:.4f}")
open(WD if False else "/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/23_亚种归属_3KRG_A12/proj/zn65_assignment.txt","w").write(
  f"ZN65 nearest subpop (30-NN vote): {c.most_common(1)[0]}\nNearest centroid: {cd[0]}\nPC1-5: {list(np.round(z,3))}\n")
print("\nVERDICT: ZN65 nearest subpopulation =", c.most_common(1)[0][0])
PY
echo "FINISH DONE $(date)"; touch "$WD/proj/FINISH_DONE"
