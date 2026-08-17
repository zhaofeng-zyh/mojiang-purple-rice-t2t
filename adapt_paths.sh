#!/usr/bin/env bash
# adapt_paths.sh — rewrite the authors' absolute paths to your own before running any pipeline step.
#
# WHY THIS EXISTS
#   The scripts under pipeline/steps/ were written to run on the authors' machine and therefore
#   contain absolute paths (/mnt/data2/...). They are published as a faithful record of how the
#   analyses were actually executed, so the paths are NOT rewritten in place upstream.
#   This helper lets you point them at your own directory layout in one step.
#
# WHY LITERAL SUBSTITUTION (and not a ${ROOT} variable)
#   Several scripts embed Python via quoted here-documents (<<'PY'). Inside a quoted here-doc the
#   shell performs NO parameter expansion, so replacing the path with ${SOMETHING} would leave the
#   literal text "${SOMETHING}" for Python to open, silently breaking those steps. Verified in this
#   repository: 10 quoted here-docs, 6 of which also contain absolute paths. Literal-to-literal
#   substitution is immune to that problem.
#
# USAGE
#   bash adapt_paths.sh --check                 # report what would change, touch nothing
#   bash adapt_paths.sh /data/my_zn65_project   # rewrite in place (a timestamped backup is made)
#   OLD_ROOT=... bash adapt_paths.sh <new_root> # override the path being replaced
#
# AFTER RUNNING
#   Your directory must mirror the layout the scripts expect, e.g.
#     <root>/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa
#     <root>/07_分析_Os02g基因鉴定/...
#   Data files themselves are not in this repository; see the Data availability statement in the
#   manuscript for the NGDC/GSA accessions.

set -euo pipefail

OLD_ROOT="${OLD_ROOT:-/mnt/data2/墨江紫米研究}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

target_files() {
  # NOTE: this script must exclude ITSELF. Its own header documents OLD_ROOT literally, so a
  # self-rewrite would clobber the default on line ~28 and make a second run unable to find the
  # original path. (Caught by running --check during development.)
  grep -rl --binary-files=without-match -F "$OLD_ROOT" \
       --include='*.sh' --include='*.py' --include='*.R' --include='*.md' . 2>/dev/null \
    | grep -vx './adapt_paths.sh' | sort
}

if [ "${1:-}" = "--check" ] || [ $# -eq 0 ]; then
  n_files=$(target_files | wc -l | tr -d ' ')
  n_lines=$(grep -r --binary-files=without-match -F -c "$OLD_ROOT" \
            --include='*.sh' --include='*.py' --include='*.R' --include='*.md' . 2>/dev/null \
            | awk -F: '{s+=$NF} END{print s+0}')
  echo "OLD_ROOT : $OLD_ROOT"
  echo "affected : $n_files files, $n_lines lines"
  target_files | sed 's/^/  /'
  [ $# -eq 0 ] && { echo; echo "Nothing was changed. Run: bash adapt_paths.sh <your_root>"; }
  exit 0
fi

NEW_ROOT="$1"
case "$NEW_ROOT" in
  /*) ;;
  *) echo "error: <new_root> must be an absolute path (got: $NEW_ROOT)" >&2; exit 2 ;;
esac
[ -d "$NEW_ROOT" ] || echo "warning: $NEW_ROOT does not exist yet; continuing anyway." >&2
case "$NEW_ROOT" in
  */) echo "error: drop the trailing slash from <new_root>" >&2; exit 2 ;;
esac

mapfile -t FILES < <(target_files)
if [ "${#FILES[@]}" -eq 0 ]; then
  echo "No file contains $OLD_ROOT — nothing to do."; exit 0
fi

STAMP=$(date +%Y%m%d_%H%M%S)
BACKUP="$HERE/.adapt_paths_backup_$STAMP"
mkdir -p "$BACKUP"
for f in "${FILES[@]}"; do
  mkdir -p "$BACKUP/$(dirname "$f")"
  cp -p "$f" "$BACKUP/$f"
done
echo "backup  : $BACKUP  (${#FILES[@]} files)"

# Literal substitution. Uses a NUL-safe python pass rather than sed so that the replacement text is
# never interpreted (paths may contain characters that are special to sed's regex/replacement).
python3 - "$OLD_ROOT" "$NEW_ROOT" "${FILES[@]}" <<'PY'
import sys
old, new, *files = sys.argv[1:]
total = 0
for path in files:
    with open(path, encoding='utf-8') as fh:
        text = fh.read()
    hits = text.count(old)
    if not hits:
        continue
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(text.replace(old, new))
    total += hits
    print("  %-58s %d" % (path, hits))
print("replaced: %d occurrences" % total)
PY

echo
echo "Syntax check:"
rc=0
for f in "${FILES[@]}"; do
  case "$f" in
    *.sh) bash -n "$f" 2>/dev/null || { echo "  FAILED (bash -n): $f" >&2; rc=1; } ;;
    *.py) python3 -m py_compile "$f" 2>/dev/null || { echo "  FAILED (py_compile): $f" >&2; rc=1; } ;;
  esac
done
find . -name '__pycache__' -prune -exec rm -rf {} + 2>/dev/null || true
if [ "$rc" -eq 0 ]; then
  echo "  all rewritten .sh/.py files still parse."
else
  echo "  SOME FILES FAILED TO PARSE — restore with: cp -a $BACKUP/. $HERE/" >&2
  exit 1
fi
echo
echo "Done. To undo: cp -a $BACKUP/. $HERE/"
