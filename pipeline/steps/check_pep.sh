#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
P=/mnt/data2/墨江紫米研究
PEP=$P/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.pep
NIP=$P/07_分析_Os02g基因鉴定/02_参考序列/nip_pep.fa
echo "=== ZN65 pep: first seq alphabet (should have many non-ACGTN = protein) ==="
seqkit head -n1 $PEP | seqkit seq -s | head -c 120; echo
echo "ZN65 non-ACGTN fraction:"; seqkit head -n 200 $PEP | seqkit seq -s | tr -d 'ACGTNacgtn\n' | wc -c | xargs echo "  non-ACGTN chars in 200 seqs:"
echo "=== NIP pep first seq ==="
seqkit head -n1 $NIP | seqkit seq -s | head -c 120; echo
echo "NIP seq count:"; seqkit stats $NIP 2>/dev/null | tail -1
