#!/bin/bash
W=/mnt/data2/墨江紫米研究/14_转录组_重比对ZN65T2T
echo "=== output tree ==="
ls -la "$W"
echo "--- counts/ ---"; ls -la "$W/counts/"
echo "--- bam/ ---"; ls "$W/bam/"
echo "--- matrix rows ---"; wc -l "$W/counts/ZN65_gene_count_matrix.tsv"
echo "--- README ---"; ls -la "$W/README_重比对说明.md"
