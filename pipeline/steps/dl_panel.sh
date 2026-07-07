WD="/Volumes/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12"
cd "$WD/genomes"
PAIRS="GCA_055776245.1:CempoIreng_black GCA_001952365.2:N22_aus GCA_000817225.1:Orufipogon_wild"
for pair in $PAIRS; do
  acc=${pair%%:*}; nm=${pair##*:}
  ok=0
  for try in 1 2 3 4 5; do
    echo "[$(date +%H:%M:%S)] $acc ($nm) try $try"
    rm -f "$acc.zip"
    if /tmp/datasets download genome accession $acc --include genome --no-progressbar --filename "$acc.zip" 2>>"$WD/logs/dl_$acc.log"; then
      if unzip -t "$acc.zip" >/dev/null 2>&1; then
        unzip -o -q "$acc.zip" -d "dir_$acc"
        fa=$(find "dir_$acc" -name '*.fna' | head -1)
        if [ -n "$fa" ]; then cp "$fa" "$nm.fna"; echo "[$(date +%H:%M:%S)] $nm OK $(du -h "$nm.fna"|cut -f1)"; ok=1; rm -f "$acc.zip"; break; fi
      fi
    fi
    echo "[$(date +%H:%M:%S)] $acc try $try failed, retrying"; sleep 5
  done
  [ $ok -eq 0 ] && echo "[$(date +%H:%M:%S)] !!! $acc FAILED after retries"
done
echo "ALL_DONE $(date +%H:%M:%S)"
ls -la "$WD/genomes/"*.fna 2>/dev/null
