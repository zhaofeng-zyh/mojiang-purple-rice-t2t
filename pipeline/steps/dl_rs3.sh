#!/bin/bash
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/10_多基因组Kala4
mkdir -p $W; cd $W
B=http://rice.hzau.edu.cn/rice_rs3/download_ext
for f in MH63RS3.fasta.gz MH63RS3_repeat.gff3.gz ZS97RS3.fasta.gz ZS97RS3_repeat.gff3.gz; do
  echo "downloading $f ..."
  wget -q --show-progress -O $f "$B/$f" 2>&1 | tail -1
  echo "  $f -> $(du -h $f 2>/dev/null | cut -f1)"
done
echo "DL_RS3_DONE"
ls -la *.gz
