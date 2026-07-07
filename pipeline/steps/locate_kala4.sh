#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
P=/mnt/data2/墨江紫米研究
W=$P/07_分析_Os02g基因鉴定/10_多基因组Kala4; cd $W
G=$P/01_基因组_组装与注释_Genome/02.annotation/02.gene
# ZN65 OsB2 protein
seqkit grep -r -n -p "ZN654G2687" $G/ZN65.longest.pep 2>/dev/null | seqkit head -n1 > osb2_zn65.pep.fa
echo "OsB2 query: $(seqkit fx2tab -nl osb2_zn65.pep.fa | cut -f2) aa"
for g in MH63RS3 ZS97RS3; do
  [ -f ${g}.ndb ] || makeblastdb -in ${g}.fasta -dbtype nucl -out $g >/dev/null 2>&1
  echo "=== OsB2/Kala4 in $g (tblastn) ==="
  tblastn -query osb2_zn65.pep.fa -db $g -outfmt "6 sseqid sstart send pident length evalue" -max_target_seqs 3 -num_threads 12 2>/dev/null | sort -k5 -rn | head -3
done
