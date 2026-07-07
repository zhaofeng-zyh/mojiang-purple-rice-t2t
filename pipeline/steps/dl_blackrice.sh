#!/bin/bash
W=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/10_多基因组Kala4; cd $W
# Cempo Ireng (black rice) GCA_055776245.1 — test NCBI vs ENA speed
URL_NCBI="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/055/776/245/GCA_055776245.1_ASM5577624v1/GCA_055776245.1_ASM5577624v1_genomic.fna.gz"
echo "trying NCBI..."
timeout 25 curl -sS -o cempoireng.fna.gz.part -w "ncbi_speed=%{speed_download}B/s got=%{size_download}\n" "$URL_NCBI" 2>&1 | tail -1
sz=$(stat -c%s cempoireng.fna.gz.part 2>/dev/null || echo 0)
echo "ncbi got $sz bytes in 25s"
# if slow (<1MB in 25s), the wrapper just records; do full download in background
rm -f cempoireng.fna.gz.part
echo "starting full NCBI download (background-safe)..."
wget -q -O cempoireng.fna.gz "$URL_NCBI" && echo "BLACKRICE_DL_DONE size=$(du -h cempoireng.fna.gz|cut -f1)" || echo "BLACKRICE_DL_FAIL"
