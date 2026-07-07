#!/bin/bash
base="http://rice.hzau.edu.cn/rice_rs3/"
echo "=== homepage download-ish links ==="
curl -sS -m 15 "$base" 2>/dev/null | grep -oiE 'href="[^"]*"' | grep -iE 'download|genome|fasta|\.gz|data|seq' | sort -u | head -30
echo ""
echo "=== candidate pages ==="
for p in download.html download.php download Download/ data.html data/ help/download.html download/index.html; do
  code=$(curl -sS -m 10 -o /dev/null -w "%{http_code}" "$base$p" 2>/dev/null)
  echo "$p -> $code"
done
echo ""
echo "=== all hrefs on homepage (first 40) ==="
curl -sS -m 15 "$base" 2>/dev/null | grep -oiE 'href="[^"]*"' | sort -u | head -40
