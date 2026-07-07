#!/bin/bash
echo "SERVER TIME: $(date '+%Y-%m-%d %H:%M:%S')"
echo "--- active jobs ---"
ps -eo etime,pcpu,comm --sort=-pcpu | grep -iE "LTR_FINDER|ltr_finder|LTR_retr|gt|diamond|orthofinder|famsa|fasttree|mcl|RepeatMask|rmblast" | grep -v grep | head -6
echo "--- LAI-fix (start 17:09) ---"; tail -2 ~/fix_lai.log 2>/dev/null
echo "--- OrthoFinder rerun ---"; tail -2 ~/rerun_of.log 2>/dev/null
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/17_batch2_质量与注释
echo "--- OF output progress ---"; ls -d $W/orthofinder/prot/OrthoFinder/Results_*/Orthogroups 2>/dev/null && echo "orthogroups produced" || echo "OF still in search/cluster"
echo "--- LAI2 files ---"; ls -lat $W/ltr2/ 2>/dev/null | head -4
