#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
M=/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript
L=/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/17_batch2_质量与注释/ltr2
echo "WAIT START: $(date '+%Y-%m-%d %H:%M:%S')"
# wait for tectonic
until [ -x ~/miniconda3/envs/tex/bin/tectonic ]; do sleep 15; done
echo "tectonic ready: $(date '+%H:%M:%S')"
conda activate tex; cd $M
tectonic manuscript.tex > compile.log 2>&1 && echo "COMPILE OK ($(date '+%H:%M:%S'))" || echo "COMPILE FAIL"
grep -ic undefined manuscript.log 2>/dev/null | xargs echo "undefined refs:"
ls -la manuscript.pdf
# wait for LAI-fix
echo "--- waiting LAI-fix ---"
until grep -q LAI_FIX_DONE ~/fix_lai.log 2>/dev/null; do sleep 20; done
echo "=== DEFINITIVE LAI (LTRharvest+LTR_FINDER) ==="
head -2 $L/genome.fa.mod.out.LAI 2>/dev/null
echo "FINALIZE DONE: $(date '+%Y-%m-%d %H:%M:%S')"
echo "FINALIZE_BATCH2_DONE"
