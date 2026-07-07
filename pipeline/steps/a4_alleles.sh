#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
P=/mnt/data2/墨江紫米研究
G=$P/01_基因组_组装与注释_Genome/02.annotation/02.gene
NIP=$P/07_分析_Os02g基因鉴定/02_参考序列/nip_pep.fa
echo "=== ZN65 vs Nipponbare protein length at pigmentation loci ==="
printf "%-14s %-16s %-8s %-8s %s\n" "locus" "ZN65_gene" "ZN65aa" "NIPaa" "note"
declare -A zng=( [Kala3_OsC1]=ZN656G0716 [Kala1_OsDFR]=ZN651G2772 [Kala4_OsB2]=ZN654G2687 [Rc]=ZN657G0823 )
declare -A nipt=( [Kala3_OsC1]=Os06t0205100 [Kala1_OsDFR]=Os01t0633500 [Kala4_OsB2]=Os04t0557500 [Rc]=Os07t0211500 )
for k in Kala3_OsC1 Kala1_OsDFR Kala4_OsB2 Rc; do
  zl=$(seqkit grep -r -n -p "${zng[$k]}" $G/ZN65.longest.pep 2>/dev/null | seqkit fx2tab -l 2>/dev/null | awk '{print $NF}' | head -1)
  nl=$(seqkit grep -r -n -p "${nipt[$k]}" $NIP 2>/dev/null | seqkit fx2tab -l 2>/dev/null | sort -k2 -rn | awk '{print $NF}' | head -1)
  printf "%-14s %-16s %-8s %-8s\n" "$k" "${zng[$k]}" "${zl:-NA}" "${nl:-NA}"
done
echo ""
echo "=== Rc focus: ZN65 Rc vs Nipponbare Rc (Nipponbare=white=rc, expect truncated) ==="
seqkit grep -r -n -p "ZN657G0823" $G/ZN65.longest.pep 2>/dev/null | seqkit seq -w0 | head -2 | tail -1 | awk '{print "ZN65 Rc protein len="length($0)}'
seqkit grep -r -n -p "Os07t0211500" $NIP 2>/dev/null | seqkit seq -w0 | head -2 | tail -1 | awk '{print "NIP Rc protein len="length($0)}'
echo "A4_DONE"
