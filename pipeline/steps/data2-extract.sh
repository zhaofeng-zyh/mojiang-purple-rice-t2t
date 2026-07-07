#!/bin/bash
set -u
ROOT=/mnt/data2/MacBookPro
W=~/data2-extract; mkdir -p "$W"
LOG="$W/extract.log"; MAN="$W/extracted_manifest.tsv"; FAILF="$W/extract_failed.txt"
: > "$LOG"; : > "$MAN"; : > "$FAILF"
log(){ echo "[$(date +%H:%M:%S)] $*" >> "$LOG"; }
# 列出真压缩包(排除._*), 按大小升序(小包先)
find "$ROOT" -type f \( -iname '*.zip' -o -iname '*.tar.gz' -o -iname '*.tgz' -o -iname '*.tar' -o -iname '*.gz' -o -iname '*.rar' -o -iname '*.7z' -o -iname '*.bz2' -o -iname '*.tar.bz2' \) ! -name '._*' -printf '%s\t%p\n' 2>/dev/null | sort -n | cut -f2- > "$W/arch_list.txt"
total=$(wc -l < "$W/arch_list.txt")
log "共 $total 个压缩包待解压"
ok=0; fail=0; skip=0; i=0
while IFS= read -r a; do
  i=$((i+1)); base=$(basename "$a"); dir=$(dirname "$a"); low="${base,,}"
  case "$low" in *.part[2-9].rar|*.part[1-9][0-9].rar) log "SKIP分卷续:$base"; skip=$((skip+1)); continue;; esac
  case "$low" in
    *.tar.gz) stem="${base%.tar.gz}";; *.tar.bz2) stem="${base%.tar.bz2}";;
    *.tgz) stem="${base%.tgz}";; *.tar) stem="${base%.tar}";; *) stem="${base%.*}";;
  esac
  ed="$dir/${stem}_解压"; mkdir -p "$ed"; rc=1
  case "$low" in
    *.zip)          unzip -o -q "$a" -d "$ed" >>"$LOG" 2>&1; rc=$?;;
    *.tar.gz|*.tgz) tar xzf "$a" -C "$ed" >>"$LOG" 2>&1; rc=$?;;
    *.tar.bz2)      tar xjf "$a" -C "$ed" >>"$LOG" 2>&1; rc=$?;;
    *.tar)          tar xf "$a" -C "$ed" >>"$LOG" 2>&1; rc=$?;;
    *.gz)           gunzip -c "$a" > "$ed/$stem" 2>>"$LOG"; rc=$?;;
    *.bz2)          bunzip2 -c "$a" > "$ed/$stem" 2>>"$LOG"; rc=$?;;
    *.rar)          unrar x -o+ -y "$a" "$ed/" >>"$LOG" 2>&1; rc=$?;;
    *.7z)           7z x -y -o"$ed" "$a" >>"$LOG" 2>&1; rc=$?;;
  esac
  if [ "$rc" -eq 0 ]; then ok=$((ok+1)); printf '%s\t%s\n' "$a" "$ed" >> "$MAN"; else fail=$((fail+1)); echo "$a" >> "$FAILF"; log "FAIL($rc):$base"; fi
  [ $((i % 10)) -eq 0 ] && log "进度 $i/$total (ok=$ok fail=$fail)"
done < "$W/arch_list.txt"
log "完成: ok=$ok fail=$fail skip=$skip"
echo "OK=$ok FAIL=$fail SKIP=$skip" > "$W/EXTRACT_DONE"
