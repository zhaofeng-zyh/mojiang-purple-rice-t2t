#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qc
P=/mnt/data2/墨江紫米研究
ASM=$P/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa
NGS=$P/00_原始测序数据_RawSequencing/基因组DNA_Genome
PEP=$P/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.pep
NIP=$P/07_分析_Os02g基因鉴定/02_参考序列/nip_pep.fa
W=$P/07_分析_Os02g基因鉴定/17_batch2_质量与注释; mkdir -p $W/{merqury,ltr,orthofinder,logs}; cd $W
echo "[$(date +%T)] A2 Merqury (meryl k=21, 12 threads to leave room for ONT)"
cd $W/merqury
meryl k=21 count output zn65_dnbseq.meryl $NGS/ZN65-1-ngs_R1.fq.gz $NGS/ZN65-1-ngs_R2.fq.gz threads=12 > ../logs/meryl.log 2>&1
merqury.sh zn65_dnbseq.meryl $ASM zn65_merqury > ../logs/merqury.log 2>&1 || echo merqury_warn
echo "  QV:"; cat zn65_merqury.qv 2>/dev/null
cd $W
echo "[$(date +%T)] A8 LTR_retriever"
cd $W/ltr; ln -sf $ASM genome.fa
gt suffixerator -db genome.fa -indexname genome -tis -suf -lcp -des -ssp -sds -dna > ../logs/gt.log 2>&1
gt ltrharvest -index genome -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seqids yes > genome.harvest.scn 2>../logs/harvest.log
LTR_retriever -genome genome.fa -inharvest genome.harvest.scn -threads 12 > ../logs/ltr.log 2>&1 || echo ltr_warn
echo "  LAI:"; cat genome.fa.out.LAI 2>/dev/null | head
cd $W
echo "[$(date +%T)] A11 OrthoFinder (ZN65 vs Nipponbare)"
cd $W/orthofinder; mkdir -p prot
seqkit seq -i $PEP > prot/ZN65.fa 2>/dev/null
seqkit seq -i $NIP > prot/Nipponbare.fa 2>/dev/null
orthofinder -f prot -t 12 -og > ../logs/of.log 2>&1 || echo of_warn
cd $W; echo "BATCH2_QC_DONE"
