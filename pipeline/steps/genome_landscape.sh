#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
P=/mnt/data2/墨江紫米研究
A=$P/01_基因组_组装与注释_Genome
W=$P/07_分析_Os02g基因鉴定/11_基因组景观; mkdir -p $W; cd $W
ASM=$A/01.assembly/ZN65.T2T.fa
GFF=$A/02.annotation/02.gene/ZN65.gff3
REP=$A/02.annotation/01.repeat/all.gff
samtools faidx $ASM 2>/dev/null
cut -f1,2 $ASM.fai | sort -V > genome.txt
# 500kb windows
bedtools makewindows -g genome.txt -w 500000 > win.bed 2>/dev/null
# gene density (genes per window)
awk -F'\t' '$3=="gene"{print $1"\t"$4"\t"$5}' $GFF | sort -k1,1 -k2,2n > genes.bed
bedtools coverage -a win.bed -b genes.bed -counts 2>/dev/null > gene_density.txt
# TE coverage fraction per window
awk -F'\t' '$3=="Transposon" || /Transposon/{print $1"\t"$4"\t"$5}' $REP 2>/dev/null | sort -k1,1 -k2,2n > te.bed
bedtools coverage -a win.bed -b te.bed 2>/dev/null > te_cov.txt
# GC per window
bedtools nuc -fi $ASM -bed win.bed 2>/dev/null | awk 'NR>1{print $1"\t"$2"\t"$3"\t"$5}' > gc.txt
echo "windows: $(wc -l < win.bed)  genes: $(wc -l < genes.bed)  TE_feats: $(wc -l < te.bed)"
echo "LANDSCAPE_DATA_DONE"
