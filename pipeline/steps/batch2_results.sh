#!/bin/bash
echo "REPORT TIME: $(date '+%Y-%m-%d %H:%M:%S')"
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/17_batch2_质量与注释
echo "=== A8 LTR Assembly Index (LAI) ==="
cat $W/ltr/genome.fa.mod.out.LAI 2>/dev/null | head -5
echo ""
echo "=== A8 intact LTR-RT superfamily counts ==="
cat $W/ltr/genome.fa.mod.out.superfam.size.list 2>/dev/null | head
echo ""
echo "=== A11 OrthoFinder summary ==="
of=$(ls -d $W/orthofinder/prot/OrthoFinder/Results_* 2>/dev/null | head -1)
echo "results dir: $of"
grep -E "Number of (species|genes|orthogroups|ortholog)" $of/Comparative_Genomics_Statistics/Statistics_Overall.tsv 2>/dev/null | head
echo "--- genes in orthogroups per species ---"
head -20 $of/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv 2>/dev/null | head -8
