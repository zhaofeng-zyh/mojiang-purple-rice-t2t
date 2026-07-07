#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
WORK=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/10_多基因组Kala4
mkdir -p $WORK; cd $WORK
B=http://ftp.ensemblgenomes.org/pub/plants/release-58/fasta
# MH63 (indica) and ZS97 (indica) and N22 (aus, sometimes pigmented) and IR64
for nm in oryza_sativa_mh63 oryza_sativa_zs97 oryza_sativa_n22 oryza_sativa_ir64; do
  short=$(echo $nm | sed 's/oryza_sativa_//')
  f=$(ls ${short}.fa.gz 2>/dev/null)
  if [ ! -s ${short}.fa.gz ] && [ ! -s ${short}.fa ]; then
    url=$(curl -sS -m 30 "$B/$nm/dna/" 2>/dev/null | grep -oE "Oryza_sativa_[^\"]*dna.toplevel.fa.gz" | head -1)
    if [ -n "$url" ]; then
      echo "downloading $short: $url"
      wget -q -O ${short}.fa.gz "$B/$nm/dna/$url" && echo "  $short OK ($(du -h ${short}.fa.gz|cut -f1))" || echo "  $short FAIL"
    else echo "$short: no url found"; fi
  fi
done
echo "=== downloaded genomes ==="
ls -la *.fa.gz 2>/dev/null
