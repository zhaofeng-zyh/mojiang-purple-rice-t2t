#!/bin/bash
# ⛔ SUPERSEDED —— 本脚本按**已被推翻的口径**输出，入口硬失败（2026-08-01 加装）。
# 旧口径：ZN65 特异序列 **73.3 Mb**。现行：**73.2 Mb**
#   = SyRI NOTAL 70,339,400 + INS 2,890,667 = **73,230,067 bp**
#   （订正见 `_口径统一_20260726/口径订正报告.md:22`）。
# 本脚本第 9、16 行的 echo 标签仍写 73.3 Mb，保留仅作留痕。
# ⚠ 为什么此前没被发现：当初给 14 个脚本加闸门针对的是 OsB2 口径，
#   而 73.3 Mb 这个数**没有任何审计模式覆盖** —— round58 加模式后才浮出（「陷阱四」第二例）。
if [ "${ALLOW_SUPERSEDED:-}" != "1" ]; then
  echo "⛔ 本脚本的 73.3 Mb 标签已被推翻（现行 73.2 Mb）。确需留痕重现请设 ALLOW_SUPERSEDED=1。" >&2
  exit 1
fi
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
