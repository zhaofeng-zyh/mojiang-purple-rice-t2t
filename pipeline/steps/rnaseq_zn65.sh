#!/bin/bash
set -e
source ~/miniconda3/etc/profile.d/conda.sh; conda activate rnaseq
P=/mnt/data2/墨江紫米研究
ASM=$P/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa
GTF=$P/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.gtf
RNA=$P/00_原始测序数据_RawSequencing/转录组RNA_Transcriptome
WORK=$P/14_转录组_重比对ZN65T2T          # NEW dir, does not overwrite originals
mkdir -p $WORK/{index,bam,counts,logs}
cd $WORK
echo "[$(date +%T)] extracting splice sites..."
hisat2_extract_splice_sites.py $GTF > index/ZN65.ss.txt
hisat2_extract_exons.py $GTF > index/ZN65.exon.txt
echo "[$(date +%T)] building HISAT2 index..."
hisat2-build -p 20 --ss index/ZN65.ss.txt --exon index/ZN65.exon.txt $ASM index/ZN65 > logs/build.log 2>&1
echo "[$(date +%T)] index done; aligning samples..."
for s in SD1 SD2 SD3 TG1 TG3; do
  echo "[$(date +%T)]  aligning $s ..."
  hisat2 -p 20 --dta -x index/ZN65 \
    -1 $RNA/$s/${s}_1.fq.gz -2 $RNA/$s/${s}_2.fq.gz \
    --new-summary --summary-file logs/${s}.summary 2> logs/${s}.hisat2.log \
    | samtools sort -@ 8 -o bam/${s}.sorted.bam -
  samtools index bam/${s}.sorted.bam
done
echo "[$(date +%T)] featureCounts..."
featureCounts -p --countReadPairs -T 20 -t exon -g gene_id -a $GTF \
  -o counts/ZN65_gene_counts.txt bam/SD1.sorted.bam bam/SD2.sorted.bam bam/SD3.sorted.bam bam/TG1.sorted.bam bam/TG3.sorted.bam > logs/featurecounts.log 2>&1
# tidy count matrix (gene + 5 samples) + mapping summary
python3 - <<PY
import re
rows=[]
with open("counts/ZN65_gene_counts.txt") as f:
    for ln in f:
        if ln.startswith("#") or ln.startswith("Geneid"): continue
        p=ln.rstrip().split("\t"); rows.append((p[0], p[6:]))
hdr=["Geneid","SD1","SD2","SD3","TG1","TG3"]
with open("counts/ZN65_gene_count_matrix.tsv","w") as o:
    o.write("\t".join(hdr)+"\n")
    for g,v in rows: o.write(g+"\t"+"\t".join(v)+"\n")
print("count matrix genes:", len(rows))
PY
echo "=== overall alignment rate per sample ==="
for s in SD1 SD2 SD3 TG1 TG3; do
  rate=$(grep -i "Overall alignment rate" logs/${s}.hisat2.log | awk "{print \$NF}")
  echo "  $s : ${rate:-NA}"
done
echo "RNASEQ_ZN65_DONE"
