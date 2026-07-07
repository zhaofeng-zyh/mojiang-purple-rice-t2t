#!/bin/bash
echo "=== reference genomes on server ==="
find /mnt/data2/墨江紫米研究 -path "*99_待审查*" -prune -o \( -iname "*nip*.fa" -o -iname "*nippon*" -o -iname "IRGSP*" -o -iname "*.fasta" \) -print 2>/dev/null | grep -iE "fa$|fasta$" | head
echo "=== existing SV analysis dirs ==="
find /mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定 -maxdepth 1 -type d | grep -iE "SV|花青素位点|比较" | head
echo "=== any existing minimap2/syri outputs ==="
find /mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定 -iname "*.delta" -o -iname "syri.out" -o -iname "*nip*.paf" 2>/dev/null | head
