WD="/Volumes/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12"
cd "$WD/genomes"
ITEMS="CempoIreng_black:https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/055/776/245/GCA_055776245.1_ASM5577624v1/GCA_055776245.1_ASM5577624v1_genomic.fna.gz N22_aus:https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/952/365/GCA_001952365.2_ASM195236v2/GCA_001952365.2_ASM195236v2_genomic.fna.gz Orufipogon_wild:https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/817/225/GCA_000817225.1_OR_W1943/GCA_000817225.1_OR_W1943_genomic.fna.gz"
for it in $ITEMS; do
  nm=${it%%:*}; url=${it#*:}
  echo "[$(date +%H:%M:%S)] $nm start"
  for try in $(seq 1 12); do
    /usr/bin/curl -sL -C - --retry 5 --retry-delay 5 --max-time 600 -o "$nm.fna.gz" "$url" && \
      /usr/bin/gzip -t "$nm.fna.gz" 2>/dev/null && { echo "[$(date +%H:%M:%S)] $nm OK $(/usr/bin/du -h "$nm.fna.gz"|cut -f1)"; break; }
    echo "[$(date +%H:%M:%S)] $nm try $try incomplete, resuming"; sleep 4
  done
done
echo "[$(date +%H:%M:%S)] decompress"
for f in CempoIreng_black N22_aus Orufipogon_wild; do
  [ -f "$f.fna.gz" ] && /usr/bin/gunzip -kf "$f.fna.gz" && echo "$f.fna $(/usr/bin/du -h "$f.fna"|cut -f1)"
done
echo "ALL_DONE $(date +%H:%M:%S)"
