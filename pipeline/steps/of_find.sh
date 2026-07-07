#!/bin/bash
of=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/17_batch2_质量与注释/orthofinder/prot/OrthoFinder/Results_Jun21
echo "=== all files produced ==="
find $of -type f | grep -iE "Orthogroups|Statistics|Counts" | head
echo "=== GeneCount header + ZN65/NIP orthogroup summary ==="
gc=$of/Orthogroups/Orthogroups.GeneCount.tsv
head -1 $gc 2>/dev/null
awk -F'\t' 'NR>1{zn=$2;nip=$3; tot++; if(zn>0&&nip>0)sh++; else if(zn>0)zo++; else no++; znsum+=zn; nipsum+=nip} END{printf "total OG=%d  shared=%d  ZN65-specific=%d  NIP-specific=%d\nZN65 genes in OG=%d  NIP genes in OG=%d\n",tot,sh,zo,no,znsum,nipsum}' $gc 2>/dev/null
echo "=== unassigned (species-specific singletons) ==="
wc -l $of/Orthogroups/Orthogroups_UnassignedGenes.tsv 2>/dev/null
