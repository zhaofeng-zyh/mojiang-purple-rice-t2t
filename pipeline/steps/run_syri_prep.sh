#!/bin/bash
set -e
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
P=/mnt/data2/墨江紫米研究
ASM=$P/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa
WORK=$P/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP
mkdir -p $WORK; cd $WORK

# 1) Download Nipponbare genome if needed
if [ ! -s nip_genome.fa ]; then
  echo "downloading Nipponbare IRGSP-1.0 ..."
  wget -q -O nip.fa.gz http://ftp.ensemblgenomes.org/pub/plants/release-58/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz
  gunzip -f nip.fa.gz
  mv Oryza_sativa.IRGSP-1.0.dna.toplevel.fa nip_raw.fa 2>/dev/null || mv nip.fa nip_raw.fa 2>/dev/null
  ls -la nip_raw.fa
fi

# 2) Extract 12 chromosomes, rename to Chr1..Chr12 (Ensembl names them 1..12)
if [ ! -s nip_genome.fa ]; then
  > nip_genome.fa
  for c in 1 2 3 4 5 6 7 8 9 10 11 12; do
    samtools faidx nip_raw.fa $c 2>/dev/null | sed "1s/^>.*/>Chr$c/" >> nip_genome.fa
  done
  samtools faidx nip_genome.fa
fi
echo "=== Nipponbare chromosomes ==="
cut -f1,2 nip_genome.fa.fai

# 3) Make a ZN65 copy with matching contig set (12 chr already)
echo "=== ready: ZN65 vs NIP chromosome sizes ==="
paste <(cut -f1,2 $ASM.fai) <(cut -f1,2 nip_genome.fa.fai)
