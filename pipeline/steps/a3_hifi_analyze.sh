#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/15_长读验证_A3; cd $W
B=bam/hifi.sorted.bam
echo "=== HiFi mean depth at key loci ==="
for reg in "Chr4:27995000-28030000:OsB2/Kala4_locus" "Chr4:27999746-28012241:OsB2_intronic_Gypsy" "Chr4:28021000-28026000:OsB2_promoter_TE_cluster" "Chr6:14000000-16000000:Chr6_inversion_region"; do
  r=$(echo $reg|cut -d: -f1-2); lab=$(echo $reg|cut -d: -f3)
  d=$(samtools depth -r "$r" $B | awk '{s+=$3;n++} END{if(n>0)printf "%.1f", s/n; else print "0"}')
  mn=$(samtools depth -r "$r" $B | awk 'BEGIN{m=1e9}{if($3<m)m=$3}END{print m+0}')
  echo "  $lab ($r): mean ${d}x, min ${mn}x"
done
echo ""
echo "=== HiFi reads single-handedly SPANNING each element (start before & end after) ==="
# OsB2 promoter TE cluster 28,021,234-28,025,736 (~4.5kb)
n1=$(samtools view $B Chr4:28021234-28025736 | awk '$4<28021234 && ($4+length($10))>28025736' | wc -l)
echo "  promoter TE cluster (28,021,234-28,025,736): $n1 HiFi reads fully span"
# OsB2 intronic Gypsy 27,999,746-28,012,241 (~12.5kb) - HiFi ~19kb can span
n2=$(samtools view $B Chr4:27999746-28012241 | awk '$4<27999746 && ($4+length($10))>28012241' | wc -l)
echo "  intronic 12.5kb Gypsy (27,999,746-28,012,241): $n2 HiFi reads fully span"
echo ""
echo "=== softclip check: well-mapped (non-softclipped) reads across OsB2 promoter ==="
tot=$(samtools view -c -F0x904 $B Chr4:28021000-28026000)
echo "  primary aligned reads over promoter: $tot"
echo "A3_HIFI_DONE"
