set -e
source ~/miniconda3/etc/profile.d/conda.sh && conda activate cgsv
B="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定"
D="$B/21_WD40伙伴_A14"
PEP="/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.pep"
FUNC="/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.pep.function.xls"
MATRIX="/mnt/data2/墨江紫米研究/14_转录组_重比对ZN65T2T/counts/ZN65_gene_count_matrix.tsv"
cd "$D"

# 1) ZN65 蛋白 DB
cp "$PEP" zn65.pep
makeblastdb -in zn65.pep -dbtype prot -out zn65db >/dev/null 2>&1

# 2) AtTTG1 -> ZN65 (找直系同源)
blastp -query AtTTG1.fasta -db zn65db -evalue 1e-5 -max_target_seqs 8 \
  -outfmt '6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore' > AtTTG1_vs_ZN65.tsv
echo "=== AtTTG1 -> ZN65 top hits (qid sid pident len ... evalue bits) ==="
column -t AtTTG1_vs_ZN65.tsv | head -8

# 3) 取最佳 ZN65 命中
BEST=$(head -1 AtTTG1_vs_ZN65.tsv | cut -f2)
echo ""
echo "=== 最佳 ZN65 直系同源: $BEST ==="
seqkit grep -n -p "$BEST" zn65.pep > best_zn65.pep 2>/dev/null
echo "ZN65 候选长度:"; seqkit fx2tab -nl best_zn65.pep
echo "起始/终止残基检查:"; seqkit fx2tab -s best_zn65.pep | awk '{seq=$2; print "start="substr(seq,1,1)"  end="substr(seq,length(seq),1)"  internal_stop="(index(substr(seq,1,length(seq)-1),"*")>0?"YES":"no")}'

# 4) 反向: ZN65候选 -> 日本晴 (RBH 确认)
blastp -query best_zn65.pep -db "$B/18_通路拷贝数_A11/nipdb" -evalue 1e-5 -max_target_seqs 5 \
  -outfmt '6 qseqid sseqid pident length evalue bitscore' > ZN65best_vs_NIP.tsv
echo ""
echo "=== ZN65候选 -> 日本晴 top hits ==="
column -t ZN65best_vs_NIP.tsv | head -6

# 5) 反向: ZN65候选 -> Arabidopsis (确认回到 TTG1)
makeblastdb -in AtTTG1.fasta -dbtype prot -out atdb >/dev/null 2>&1
blastp -query best_zn65.pep -db atdb -evalue 1e-3 -outfmt '6 qseqid sseqid pident length evalue bitscore' > ZN65best_vs_At.tsv
echo ""
echo "=== ZN65候选 -> AtTTG1 (反向确认) ==="
column -t ZN65best_vs_At.tsv | head -3

# 6) 功能注释 (WD40 域?)
echo ""
echo "=== $BEST 功能注释行 ==="
grep "${BEST%.*}" "$FUNC" | head -2 | cut -c1-400

# 7) 表达 (Step14 count matrix)
GID="${BEST%.*}"
echo ""
echo "=== 表达 (Step14 count matrix) 表头 + $GID 行 ==="
head -1 "$MATRIX"
grep -P "^${GID}\t" "$MATRIX" || echo "（矩阵中未找到 $GID，尝试不带版本/前缀匹配）"
grep "$GID" "$MATRIX" | head -3
