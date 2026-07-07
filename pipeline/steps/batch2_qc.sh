#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
P=/mnt/data2/墨江紫米研究
ASM=$P/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa
NGS=$P/00_原始测序数据_RawSequencing/基因组DNA_Genome
PEP=$P/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.pep
NIP=$P/07_分析_Os02g基因鉴定/02_参考序列/nip_pep.fa
W=$P/07_分析_Os02g基因鉴定/17_batch2_质量与注释; mkdir -p $W/{merqury,ltr,orthofinder,logs}; cd $W
# wait for A3 to finish (free cores) and qc env ready
echo "[$(date +%T)] waiting for A3 + qc env..."
until grep -q A3_LONGREAD_DONE ~/a3_longread.log 2>/dev/null; do sleep 30; done
until [ -x ~/miniconda3/envs/qc/bin/merqury.sh ] || [ -x ~/miniconda3/envs/qc/bin/merqury ]; do sleep 30; done
conda activate qc

echo "[$(date +%T)] A2: Merqury QV (meryl k=21 from DNBSEQ)"
cd $W/merqury
meryl k=21 count output zn65_dnbseq.meryl $NGS/ZN65-1-ngs_R1.fq.gz $NGS/ZN65-1-ngs_R2.fq.gz threads=20 > ../logs/meryl.log 2>&1
merqury.sh zn65_dnbseq.meryl $ASM zn65_merqury > ../logs/merqury.log 2>&1 || echo "merqury warn"
echo "  QV:"; cat zn65_merqury.qv 2>/dev/null
cd $W

echo "[$(date +%T)] A8: LTR_retriever (LAI + intact LTR-RTs)"
cd $W/ltr
ln -sf $ASM genome.fa
gt suffixerator -db genome.fa -indexname genome -tis -suf -lcp -des -ssp -sds -dna > ../logs/gt_suf.log 2>&1
gt ltrharvest -index genome -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > genome.harvest.scn 2>../logs/ltrharvest.log
LTR_retriever -genome genome.fa -inharvest genome.harvest.scn -threads 20 > ../logs/ltr_retriever.log 2>&1 || echo "ltr_retriever warn"
echo "  LAI:"; cat genome.fa.out.LAI 2>/dev/null | head
cd $W

echo "[$(date +%T)] A11: OrthoFinder (ZN65 vs Nipponbare)"
cd $W/orthofinder; mkdir -p prot
# clean headers (first word) to avoid OrthoFinder issues
seqkit seq -i $PEP > prot/ZN65.fa 2>/dev/null
seqkit seq -i $NIP > prot/Nipponbare.fa 2>/dev/null
orthofinder -f prot -t 20 -og > ../logs/orthofinder.log 2>&1 || echo "orthofinder warn"
cd $W
echo "BATCH2_QC_DONE"
