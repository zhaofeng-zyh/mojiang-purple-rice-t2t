#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
P=/mnt/data2/墨江紫米研究
G=$P/01_基因组_组装与注释_Genome/02.annotation/02.gene
NIP=$P/07_分析_Os02g基因鉴定/02_参考序列/nip_pep.fa
W=$P/07_分析_Os02g基因鉴定/16_着色位点等位_A4; mkdir -p $W; cd $W
align(){ # zn_gene nip_gene label
  seqkit grep -r -n -p "$1" $G/ZN65.longest.pep | seqkit head -n1 | seqkit seq -w0 | sed "1s/.*/>ZN65_$3/" > z.fa
  # pick longest NIP isoform for this gene
  seqkit grep -r -n -p "$2" $NIP | seqkit sort -l -r 2>/dev/null | seqkit head -n1 | seqkit seq -w0 | sed "1s/.*/>NIP_$3/" > n.fa
  cat z.fa n.fa > ${3}.fa
  mafft --quiet --auto ${3}.fa > ${3}.aln 2>/dev/null
  python3 - "$3" <<PY
import sys
lab=sys.argv[1]; d={};h=None
for l in open(lab+".aln"):
    if l.startswith(">"): h=l[1:].split()[0]; d[h]=""
    else: d[h]+=l.strip()
k=list(d)
if len(k)<2: print(lab+": MISSING one seq"); sys.exit()
s1,s2=d[k[0]],d[k[1]]
ident=sum(1 for a,b in zip(s1,s2) if a==b and a!="-")
both=sum(1 for a,b in zip(s1,s2) if a!="-" and b!="-")
ntail=len(s2)-len(s2.rstrip("-")); ztail=len(s1)-len(s1.rstrip("-"))
print(f"{lab}: ZN65 {len(s1.replace('-',''))}aa vs NIP {len(s2.replace('-',''))}aa | overlap {both}aa id {100*ident/both:.1f}% | NIP missing C-term vs ZN65 = {ntail}aa ; ZN65 missing = {ztail}aa")
PY
}
align ZN657G0823 Os07g0211500 Rc
align ZN654G2687 Os04g0557500 Kala4_OsB2
align ZN656G0716 Os06g0205100 Kala3_OsC1
align ZN651G2772 Os01g0633500 Kala1_OsDFR
echo "A4_ALIGN_DONE"
