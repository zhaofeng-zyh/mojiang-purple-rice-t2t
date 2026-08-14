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
# 1) wait for download to finish
echo "waiting for black rice download..."
until grep -q BLACKRICE_DL_DONE ~/dl_blackrice.log 2>/dev/null; do sleep 20; done
echo "download done; decompressing..."
gunzip -kf cempoireng.fna.gz
samtools faidx cempoireng.fna 2>/dev/null
echo "=== Cempo Ireng (black rice) chromosomes ==="
grep ">" cempoireng.fna | head -14 | sed 's/>//' | awk '{print $1}'
# 2) locate OsB2/Kala4 in black rice via tblastn (ZN65 OsB2 protein)
[ -f cempo.ndb ] || makeblastdb -in cempoireng.fna -dbtype nucl -out cempo >/dev/null 2>&1
echo "=== OsB2/Kala4 in Cempo Ireng (tblastn) ==="
tblastn -query osb2_zn65.pep.fa -db cempo -outfmt "6 sseqid sstart send pident length evalue" -max_target_seqs 3 -num_threads 12 2>/dev/null | sort -k5 -rn | head -3
top=$(tblastn -query osb2_zn65.pep.fa -db cempo -outfmt "6 sseqid sstart send" -max_target_seqs 1 -num_threads 12 2>/dev/null | head -1)
chr=$(echo "$top" | cut -f1); s=$(echo "$top" | cut -f2); e=$(echo "$top" | cut -f3)
lo=$(( (s<e?s:e) - 25000 )); hi=$(( (s>e?s:e) + 25000 ))
echo "OsB2 locus window: $chr:$lo-$hi"
samtools faidx cempoireng.fna "$chr:$lo-$hi" 2>/dev/null | seqkit seq -w0 | sed '1s/.*/>CempoIreng_Kala4/' > CempoIreng_Kala4_locus.fasta
# 3) DECISIVE: does black rice OsB2 locus contain the ZN65 OsB2 promoter-proximal insertion?
makeblastdb -in CempoIreng_Kala4_locus.fasta -dbtype nucl -out cempoloc >/dev/null 2>&1
echo "=== ZN65 OsB2 promoter-proximal (3.5kb) vs Cempo Ireng OsB2 locus ==="
blastn -query osb2_prom_zn65.fa -db cempoloc -outfmt "6 qstart qend length pident" -evalue 1e-5 2>/dev/null | awk '{print "  qpos "$1"-"$2"  len="$3"  id="$4"%"}' | head
cov=$(blastn -query osb2_prom_zn65.fa -db cempoloc -outfmt "6 length" -evalue 1e-5 2>/dev/null | awk '{s+=$1} END{print s+0}')
echo "  -> ZN65 promoter aligned to BLACK rice: $cov bp (of 3504)"
echo "  (recall: vs white Nipponbare 984bp proximal only; vs white MH63/ZS97 = 0bp)"
# 4) whole OsB2 locus alignment ZN65 vs black rice
nucmer --maxmatch -c 100 -p znVScempo $P/07_分析_Os02g基因鉴定/07_花青素位点SV/OsB2_Kala4/ZN65_Kala4_locus.fasta CempoIreng_Kala4_locus.fasta 2>/dev/null
dnadiff -d znVScempo.delta -p znVScempo 2>/dev/null
echo "=== ZN65 vs Cempo Ireng OsB2 locus identity/aligned ==="
grep -E "AlignedBases|AvgIdentity" znVScempo.report | head -3
echo "KALA4_BLACKRICE_DONE"
