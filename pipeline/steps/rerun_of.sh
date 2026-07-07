#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qc
echo "OF START: $(date '+%Y-%m-%d %H:%M:%S')"
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/17_batch2_质量与注释/orthofinder
cd $W
rm -rf prot/OrthoFinder   # clear partial run
orthofinder -f prot -t 16 -og > of_rerun.log 2>&1 && echo OF_OK || echo OF_FAIL
echo "OF FINISH: $(date '+%Y-%m-%d %H:%M:%S')"
of=$(ls -d prot/OrthoFinder/Results_*/ 2>/dev/null | head -1)
echo "results: $of"
head -2 ${of}Comparative_Genomics_Statistics/Statistics_Overall.tsv 2>/dev/null
echo "OF_RERUN_DONE"
