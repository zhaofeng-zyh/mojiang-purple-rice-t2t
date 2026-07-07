#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP; cd $W
# Filter out base-level variants (SNP/INS/DEL/CPG/CPL) -> keep structural for plotsr
awk -F'\t' '$11!="SNP" && $11!="INS" && $11!="DEL" && $11!="CPG" && $11!="CPL"' zn65nip_syri.out > syri_structural.out
echo "structural lines: $(wc -l < syri_structural.out)"
printf "#file\tname\ttags\nnip_genome.fa\tNipponbare\tlw:1.5\nzn65_genome.fa\tZN65_purple\tlw:1.5\n" > genomes.txt
timeout 300 plotsr --sr syri_structural.out --genomes genomes.txt -o zn65_vs_nip_SV.png -H 8 -W 11 2>plotsr.err && echo PLOTSR_DONE || { echo PLOTSR_FAIL; tail -8 plotsr.err; }
ls -la zn65_vs_nip_SV.png 2>/dev/null
