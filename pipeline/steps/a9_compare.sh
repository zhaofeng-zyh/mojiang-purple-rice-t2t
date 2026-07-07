#!/bin/bash
A5=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP/zn65nip_syri.out
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/19_SV正交验证_A9
echo "=== High-level structural SV counts: asm5 vs asm10 vs asm20 ==="
printf "%-14s %-8s %-8s %-8s\n" "Type" "asm5" "asm10" "asm20"
for t in SYN INV TRANS INVTR DUP INVDP NOTAL; do
  c5=$(awk -F"\t" -v t=$t "\$11==t" $A5 2>/dev/null | wc -l)
  c10=$(awk -F"\t" -v t=$t "\$11==t" $W/asm10_syri.out 2>/dev/null | wc -l)
  c20=$(awk -F"\t" -v t=$t "\$11==t" $W/asm20_syri.out 2>/dev/null | wc -l)
  printf "%-14s %-8s %-8s %-8s\n" "$t" "$c5" "$c10" "$c20"
done
echo ""
echo "=== Assemblytics / nucmer status ==="
ls -la $W/zn65nip.delta 2>/dev/null | awk "{print \"delta: \"\$5\" bytes\"}"
ls $W/assemblytics* 2>/dev/null | head -3 || echo "no Assemblytics output"
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qc 2>/dev/null
which Assemblytics 2>/dev/null || echo "Assemblytics NOT installed"
which svim-asm 2>/dev/null && echo "svim-asm available" || echo "svim-asm not installed"
