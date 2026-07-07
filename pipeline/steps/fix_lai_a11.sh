#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qc
echo "START: $(date '+%Y-%m-%d %H:%M:%S')"
P=/mnt/data2/墨江紫米研究
ASM=$P/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa
W=$P/07_分析_Os02g基因鉴定/17_batch2_质量与注释/ltr2; mkdir -p $W; cd $W
ln -sf $ASM genome.fa
# LTRharvest
gt suffixerator -db genome.fa -indexname gt -tis -suf -lcp -des -ssp -sds -dna > gt.log 2>&1
gt ltrharvest -index gt -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seqids yes > harvest.scn 2>harvest.log
# LTR_FINDER (parallel wrapper if available, else ltr_finder)
if command -v LTR_FINDER_parallel >/dev/null 2>&1; then
  LTR_FINDER_parallel -seq genome.fa -threads 16 -harvest_out > finder.log 2>&1
  cat genome.fa.finder.combine.scn > finder.scn 2>/dev/null
fi
# combine and run LTR_retriever
if [ -s finder.scn ]; then
  cat harvest.scn finder.scn > combined.scn
  LTR_retriever -genome genome.fa -inharvest combined.scn -threads 16 > ltr.log 2>&1
else
  LTR_retriever -genome genome.fa -inharvest harvest.scn -threads 16 > ltr.log 2>&1
fi
echo "=== definitive LAI ==="; head -2 genome.fa.mod.out.LAI 2>/dev/null
echo "FINISH: $(date '+%Y-%m-%d %H:%M:%S')"
echo "LAI_FIX_DONE"
