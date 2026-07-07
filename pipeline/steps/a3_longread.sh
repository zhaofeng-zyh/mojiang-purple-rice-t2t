#!/bin/bash
set -e
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
P=/mnt/data2/墨江紫米研究
ASM=$P/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa
RAW=$P/00_原始测序数据_RawSequencing/基因组DNA_Genome
W=$P/07_分析_Os02g基因鉴定/15_长读验证_A3; mkdir -p $W/bam $W/logs; cd $W
echo "[$(date +%T)] HiFi: BAM->fastq->minimap2 map-hifi"
samtools fastq -@ 8 $RAW/ZN65-1.hifi_reads.bam 2>logs/hifi_fastq.log | \
  minimap2 -ax map-hifi -t 20 $ASM - 2>logs/hifi_mm2.log | \
  samtools sort -@ 8 -o bam/hifi.sorted.bam -; samtools index bam/hifi.sorted.bam
echo "[$(date +%T)] ONT-UL: minimap2 map-ont"
minimap2 -ax map-ont -t 20 $ASM $RAW/ZN65-1.pass.ul.fq.gz 2>logs/ont_mm2.log | \
  samtools sort -@ 8 -o bam/ont_ul.sorted.bam -; samtools index bam/ont_ul.sorted.bam
echo "[$(date +%T)] depth at key loci"
# OsB2/Kala4 locus + Chr6 inversion region (approx 13-17 Mb from plotsr)
for reg in "Chr4:27995000-28030000" "Chr6:13000000-17000000"; do
  echo "== $reg ==" >> logs/depth_summary.txt
  echo "HiFi mean depth:" >> logs/depth_summary.txt
  samtools depth -r "$reg" bam/hifi.sorted.bam | awk '{s+=$3;n++} END{printf "  %.1f (n=%d)\n", s/n, n}' >> logs/depth_summary.txt
  echo "ONT-UL mean depth:" >> logs/depth_summary.txt
  samtools depth -r "$reg" bam/ont_ul.sorted.bam | awk '{s+=$3;n++} END{printf "  %.1f (n=%d)\n", s/n, n}' >> logs/depth_summary.txt
done
# spanning reads across the OsB2 Gypsy insertion (must span 27,999,746-28,012,241, ~12.5kb)
echo "== reads spanning OsB2 Gypsy (Chr4:27,999,000-28,013,000) ==" >> logs/depth_summary.txt
samtools view bam/ont_ul.sorted.bam Chr4:27999000-28013000 | awk '$4<27999000 && ($4+length($10))>28013000' | wc -l | xargs echo "  ONT-UL reads fully spanning (single-read):" >> logs/depth_summary.txt
cat logs/depth_summary.txt
echo "A3_LONGREAD_DONE"
