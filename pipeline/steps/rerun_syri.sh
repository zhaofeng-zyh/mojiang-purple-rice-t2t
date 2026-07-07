#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
WORK=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP
cd $WORK
syri -c zn65_vs_nip.bam -r nip_genome.fa -q zn65_genome.fa -F B --prefix zn65nip_ --nc 12 2>syri.err
if [ -s zn65nip_syri.out ]; then echo SYRI_OK; else echo SYRI_ERR1; tail -5 syri.err; fi
