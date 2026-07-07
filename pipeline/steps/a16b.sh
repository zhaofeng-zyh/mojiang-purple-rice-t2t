#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
RM=/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/02.annotation/01.repeat/repeatmasker.gff
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/20_TE家族归因_A16; cd $W
# build family BED from repeatmasker.gff
awk -F"\t" '$3=="Transposon"{match($9,/Class=([^;]+)/,a); cls=a[1]; sub(/\/.*/,"",cls); fam=cls; if(a[1]~/Gypsy/)fam="LTR/Gypsy"; else if(a[1]~/Copia/)fam="LTR/Copia"; else if(a[1]~/^LTR/)fam="LTR/other"; print $1"\t"$4-1"\t"$5"\t"fam}' $RM | sort -k1,1 -k2,2n > rm_family.bed
echo "RM features: $(wc -l < rm_family.bed)"
echo "chr names check (RM vs specific):"; cut -f1 rm_family.bed | sort -u | head -3 | tr "\n" " "; echo; cut -f1 zn65_specific.bed | sort -u | head -3 | tr "\n" " "; echo
echo "=== per-family bp WITHIN ZN65-specific (73.3 Mb gained) sequence ==="
for fam in "LTR/Gypsy" "LTR/Copia" "LTR/other" "DNA" "LINE" "SINE"; do
  grep -P "\t${fam}$" rm_family.bed > f.bed
  bp=$(bedtools intersect -a zn65_specific.bed -b f.bed 2>/dev/null | sort -k1,1 -k2,2n | bedtools merge 2>/dev/null | awk "{s+=\$3-\$2} END{printf \"%.2f\", s/1e6}")
  echo "  ${fam}: ${bp:-0} Mb"
done
allrep=$(bedtools intersect -a zn65_specific.bed -b rm_family.bed 2>/dev/null | sort -k1,1 -k2,2n | bedtools merge 2>/dev/null | awk "{s+=\$3-\$2} END{printf \"%.2f\", s/1e6}")
echo "  ALL repeats within ZN65-specific: ${allrep} Mb of 73.3 Mb"
echo "A16B_DONE"
