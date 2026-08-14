#!/bin/bash
# ⛔ SUPERSEDED —— 本脚本按**已被推翻的口径**输出，入口硬失败（2026-07-29 加装）。
# 旧口径「OsB2 启动子 ~5.9/6.3 kb ZN65 特异插入、白等位缺失、存在-缺失共显性」已被推翻；
# 现行：上游 ~5.87 kb block A 复合局部片段重复，ZN65 3 拷贝 / 6 对照各 1 拷贝，
# 片段在全部 7 个基因组中都存在（差异只是拷贝数），接头缺 TSD 与末端重复、非逆转座子插入；
# LINE1-11_OS 实测 5,934 bp；三拷贝区 A1–A3 = 115,450 bp。
if [ "${ALLOW_SUPERSEDED:-}" != "1" ]; then
  echo "⛔ 本脚本按已被推翻的口径输出（详见文件头注释）。确需留痕重现请设 ALLOW_SUPERSEDED=1。" >&2
  exit 1
fi
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
P=/mnt/data2/墨江紫米研究
W=$P/07_分析_Os02g基因鉴定/10_多基因组Kala4; cd $W
ZNLOC=$P/07_分析_Os02g基因鉴定/07_花青素位点SV/OsB2_Kala4/ZN65_Kala4_locus.fasta
# extract indica OsB2 loci (60kb windows around OsB2)
seqkit faidx MH63RS3.fasta "Chr04:29605000-29665000" 2>/dev/null | seqkit seq -w0 | sed '1s/.*/>MH63_Kala4/' > MH63_Kala4_locus.fasta
seqkit faidx ZS97RS3.fasta "Chr04:29380000-29440000" 2>/dev/null | seqkit seq -w0 | sed '1s/.*/>ZS97_Kala4/' > ZS97_Kala4_locus.fasta
echo "=== large retroelements (LTR/LINE) at OsB2 locus per genome (repeat gff) ==="
echo "--- MH63 OsB2 region Chr04:29,605,000-29,665,000 ---"
zcat MH63RS3_repeat.gff3.gz 2>/dev/null | awk -F'\t' '$1=="Chr04" && $4<29665000 && $5>29605000 && ($9 ~ /LTR|Gypsy|Copia|LINE/){n++; l+=$5-$4} END{print "  LTR/LINE elements:",n+0," total bp:",l+0}'
echo "--- ZS97 OsB2 region Chr04:29,380,000-29,440,000 ---"
zcat ZS97RS3_repeat.gff3.gz 2>/dev/null | awk -F'\t' '$1=="Chr04" && $4<29440000 && $5>29380000 && ($9 ~ /LTR|Gypsy|Copia|LINE/){n++; l+=$5-$4} END{print "  LTR/LINE elements:",n+0," total bp:",l+0}'
echo "  (ZN65 OsB2 locus had: 12.5kb Gypsy RETRO2B + 6.3kb LINE1 = ~19kb retroelements)"
echo ""
echo "=== align ZN65 OsB2 locus vs MH63 / ZS97 (aligned% reveals ZN65-specific insertions) ==="
for q in MH63 ZS97; do
  nucmer --maxmatch -c 100 -p ${q}cmp $ZNLOC ${q}_Kala4_locus.fasta 2>/dev/null
  dnadiff -d ${q}cmp.delta -p ${q}cmp 2>/dev/null
  echo "--- ZN65 vs $q ---"; grep -E "AlignedBases|AvgIdentity" ${q}cmp.report | head -3
done
echo "DONE"
