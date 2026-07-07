#!/bin/bash
set -e
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
P=/mnt/data2/墨江紫米研究
ASM=$P/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa
WORK=$P/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP
cd $WORK
ln -sf $ASM zn65_genome.fa
samtools faidx zn65_genome.fa 2>/dev/null || true

echo "=== minimap2 whole-genome alignment (ref=NIP, qry=ZN65, asm5) ==="
minimap2 -ax asm5 --eqx -t 20 nip_genome.fa zn65_genome.fa 2>mm2.log | samtools sort -@ 8 -O BAM -o zn65_vs_nip.bam -
samtools index zn65_vs_nip.bam
echo "alignment done: $(du -h zn65_vs_nip.bam | cut -f1)"

echo "=== per-chromosome primary alignment strand (check for whole-chr inversions) ==="
samtools view -F 0x100 -q 10 zn65_vs_nip.bam | awk '{f=and($2,16)?"-":"+"; print $3"\t"f}' | sort | uniq -c | sort -k2 | head -30

echo "=== SyRI ==="
syri -c zn65_vs_nip.bam -r nip_genome.fa -q zn65_genome.fa -F B --prefix zn65nip_ --nc 12 2>syri.err || { echo "SYRI_ERR"; tail -20 syri.err; }
ls -la zn65nip_syri.out 2>/dev/null && echo "SYRI_OK"
