#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
echo "ASSEMBLYTICS START: $(date '+%Y-%m-%d %H:%M:%S')"
D=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/19_SV正交验证_A9; cd $W
conda activate cgsv
nucmer --maxmatch -l 100 -c 500 -t 20 -p zn65nip $D/nip_genome.fa $D/zn65_genome.fa 2>nucmer2.log
echo "delta: $(ls -la zn65nip.delta 2>/dev/null | awk "{print \$5}") bytes ($(date +%T))"
conda activate qc
Assemblytics zn65nip.delta asmtics 10000 50 10000 > asmtics.log 2>&1 && echo "ASMTICS OK" || tail -3 asmtics.log
echo "=== Assemblytics orthogonal SV summary ==="
cat asmtics.Assemblytics_structural_variants.summary 2>/dev/null | head -30
echo "ASSEMBLYTICS FINISH: $(date '+%Y-%m-%d %H:%M:%S')"
echo "ASMTICS_DONE"
