#!/bin/bash
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/17_batch2_质量与注释
echo "NOW: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=== A2 Merqury QV ==="
cat $W/merqury/zn65_merqury.qv 2>/dev/null || echo "no qv yet"
echo "=== all batch2/LTR-related processes ==="
ps -eo pid,etime,pcpu,comm,args --sort=-pcpu | grep -iE "batch2b|LTR_retriever|RepeatMask|rmblast|cd-hit|TEsorter|makeblast|blast|orthofinder|trf|HMMER|hmm" | grep -v grep | head -10
echo "=== batch2b log (with timestamps it printed) ==="
cat ~/batch2b_qc.log 2>/dev/null
echo "=== LTR log tail ==="
tail -6 $W/logs/ltr.log 2>/dev/null
echo "=== LTR dir files (mtime = growth?) ==="
ls -lat $W/ltr/ 2>/dev/null | head -8
