#!/bin/bash
G=/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/02.annotation/02.gene
echo "=== GTF exon/CDS lines (need gene_id + transcript_id) ==="
grep -P "\texon\t" "$G/ZN65.gtf" | head -2
echo "=== GTF feature types ==="
cut -f3 "$G/ZN65.gtf" | sort | uniq -c | head
echo "=== gff3 sample ==="
grep -P "\tgene\t|\tmRNA\t|\texon\t" "$G/ZN65.gff3" | head -3
echo "=== read length SD1_1 ==="
zcat /mnt/data2/墨江紫米研究/00_原始测序数据_RawSequencing/转录组RNA_Transcriptome/SD1/SD1_1.fq.gz | head -2 | tail -1 | awk '{print length($0)" bp"}'
