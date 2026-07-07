#!/bin/bash
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/17_batch2_质量与注释/orthofinder/prot/OrthoFinder
of=$(ls -dt $W/Results_*/ 2>/dev/null | head -1)
echo "results: $of"
echo "=== Overall stats ==="
grep -iE "Number of species|Number of genes$|Number of orthogroups|in orthogroups|species-specific orthogroups|G50|O50|single-copy" "$of/Comparative_Genomics_Statistics/Statistics_Overall.tsv" 2>/dev/null | head -15
echo "=== Per-species ==="
cat "$of/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv" 2>/dev/null | awk -F'\t' 'NR<=6{print}'
echo "=== LAI-fix progress ==="
tail -2 ~/fix_lai.log 2>/dev/null
ls -lat /mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/17_batch2_质量与注释/ltr2/ 2>/dev/null | head -3
