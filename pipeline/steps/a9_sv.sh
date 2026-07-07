#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
echo "A9 START: $(date '+%Y-%m-%d %H:%M:%S')"
D=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/19_SV正交验证_A9; mkdir -p $W; cd $W
REF=$D/nip_genome.fa; QRY=$D/zn65_genome.fa
for preset in asm10 asm20; do
  echo "[$(date +%T)] minimap2 $preset ..."
  minimap2 -ax $preset --eqx -t 20 $REF $QRY 2>$preset.mm2.log | samtools sort -@ 8 -O BAM - > $preset.bam 2>/dev/null
  samtools index $preset.bam
  echo "[$(date +%T)] syri $preset ..."
  syri -c $preset.bam -r $REF -q $QRY -F B --prefix ${preset}_ --nc 8 -k 2>$preset.syri.log >/dev/null || echo "$preset syri warn"
done
echo "[$(date +%T)] Assemblytics (nucmer) ..."
conda activate qc 2>/dev/null || conda activate cgsv
nucmer --maxmatch -l 100 -c 500 -t 20 -p zn65nip $REF $QRY 2>nucmer.log
if command -v Assemblytics >/dev/null 2>&1; then
  Assemblytics zn65nip.delta assemblytics 10000 50 10000 > assemblytics.log 2>&1 || echo "assemblytics warn"
fi
echo "=== SV count comparison (structural: INV/TRANS/DUP/INS/DEL) ==="
for p in asm10 asm20; do
  echo "--- $p ---"; grep -vE "SNP|CPG|CPL|HDR|TDM" ${p}_syri.out 2>/dev/null | awk "{print \$11}" | sort | uniq -c | sort -rn | head
done
echo "=== Assemblytics summary ==="
cat assemblytics.Assemblytics_structural_variants.summary 2>/dev/null | head -20
echo "A9 FINISH: $(date '+%Y-%m-%d %H:%M:%S')"
echo "A9_DONE"
