#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
echo "A16 START: $(date '+%Y-%m-%d %H:%M:%S')"
D=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP
LTR=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/17_batch2_质量与注释/ltr2/genome.fa.mod.LTR.gff3
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/20_TE家族归因_A16; mkdir -p $W; cd $W
# 1) ZN65-specific sequence (query side): NOTAL on query + INS. syri.out cols: 1RefChr 2Rs 3Re 4Rseq 5Qseq 6QryChr 7Qs 8Qe 9id 10pid 11type
awk -F"\t" '($11=="NOTAL" || $11=="INS") && $6!="-" && $7!="-" {s=$7; e=$8; if(s>e){t=s;s=e;e=t} if(e>s) print $6"\t"s-1"\t"e}' $D/zn65nip_syri.out | sort -k1,1 -k2,2n > zn65_specific.bed
echo "ZN65-specific regions: $(wc -l < zn65_specific.bed) ; total bp:"; awk '{s+=$3-$2} END{printf "  %.1f Mb\n", s/1e6}' zn65_specific.bed
# 2) intact LTR-RT annotation -> BED, with family
grep -P "\tLTR_retriever\t(repeat_region|LTR_retrotransposon)\t" $LTR 2>/dev/null | awk -F"\t" '{fam="LTR"; if($9~/Gypsy/)fam="Gypsy"; else if($9~/Copia/)fam="Copia"; print $1"\t"$4-1"\t"$5"\t"fam}' | sort -k1,1 -k2,2n > ltr_intact.bed
echo "intact LTR-RT features: $(wc -l < ltr_intact.bed)"
# 3) intersect: how much ZN65-specific seq is intact LTR-RT, by family
echo "=== ZN65-specific sequence overlapped by intact LTR-RTs (by family) ==="
for fam in Gypsy Copia LTR; do
  grep -P "\t$fam$" ltr_intact.bed > f.bed
  bp=$(bedtools intersect -a zn65_specific.bed -b f.bed 2>/dev/null | awk '{s+=$3-$2} END{printf "%.2f", s/1e6}')
  echo "  $fam : ${bp:-0} Mb of ZN65-specific seq"
done
totbp=$(bedtools intersect -a zn65_specific.bed -b ltr_intact.bed 2>/dev/null | awk '{s+=$3-$2} END{printf "%.2f", s/1e6}')
echo "  ALL intact LTR-RT : ${totbp} Mb"
echo "A16 FINISH: $(date '+%Y-%m-%d %H:%M:%S')"
echo "A16_DONE"
