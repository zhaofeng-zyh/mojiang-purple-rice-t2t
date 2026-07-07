#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
WORK=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP
cd $WORK
echo "=== SyRI output lines: $(wc -l < zn65nip_syri.out) ==="
echo "=== SV/annotation type counts (col 11) ==="
awk -F'\t' '{print $11}' zn65nip_syri.out | sort | uniq -c | sort -rn
echo ""
echo "=== Structural rearrangement summary (bp on reference) ==="
awk -F'\t' '$11 ~ /^(SYN|INV|TRANS|INVTR|DUP|INVDP|NOTAL)$/ && $2!="-" && $3!="-"{
  len=$3-$2+1; tot[$11]+=len; n[$11]++
} END{for(t in tot) printf "  %-8s %6d events  %12d bp\n", t, n[t], tot[t]}'
echo ""
echo "=== small variants ==="
for t in SNP INS DEL CPG HDR TDM; do c=$(awk -F'\t' -v T=$t '$11==T' zn65nip_syri.out | wc -l); echo "  $t: $c"; done

echo ""
echo "=== plotsr genome-wide figure ==="
printf "#file\tname\ttags\nnip_genome.fa\tNipponbare\tlw:1.5\nzn65_genome.fa\tZN65(purple)\tlw:1.5\n" > genomes.txt
plotsr --sr zn65nip_syri.out --genomes genomes.txt -o zn65_vs_nip_SV.png 2>plotsr.err && echo PLOTSR_OK || { echo PLOTSR_ERR; tail -5 plotsr.err; }
ls -la zn65_vs_nip_SV.png 2>/dev/null

echo ""
echo "=== SVs at Kala4/OsB2 locus (NIP Chr4:27,915,598-27,939,357) ==="
awk -F'\t' '$1=="Chr4" && $2!="-" && $2>27905000 && $3<27950000 && $11 ~ /^(INV|TRANS|DUP|INVTR|INVDP|NOTAL|CPG|HDR|TDM)$/{print $1":"$2"-"$3"\t"$11"\t"$6":"$7"-"$8}' zn65nip_syri.out | head -25
echo "--- OsC1 locus (NIP Chr6:5,315,178-5,316,875) ---"
awk -F'\t' '$1=="Chr6" && $2!="-" && $2>5305000 && $3<5325000 && $11 ~ /^(INV|TRANS|DUP|INVTR|INVDP|NOTAL|CPG|HDR|TDM)$/{print $1":"$2"-"$3"\t"$11"\t"$6":"$7"-"$8}' zn65nip_syri.out | head
