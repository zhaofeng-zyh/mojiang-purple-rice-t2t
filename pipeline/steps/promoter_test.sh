#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
P=/mnt/data2/墨江紫米研究
W=$P/07_分析_Os02g基因鉴定/10_多基因组Kala4; cd $W
ASM=$P/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa
# ZN65 OsB2 promoter-proximal: TSS at Chr4:28,021,147 (minus), upstream = higher coords. Take TSS..+3.5kb
seqkit faidx $ASM "Chr4:28021147-28024650" 2>/dev/null | seqkit seq -w0 | sed '1s/.*/>ZN65_OsB2_promoter_proximal/' > osb2_prom_zn65.fa
echo "ZN65 OsB2 promoter-proximal: $(seqkit fx2tab -nl osb2_prom_zn65.fa | cut -f2) bp"
# BLAST it against each genome's OsB2 LOCUS window (not whole genome, to avoid TE multi-copy noise)
NIPLOC=$P/07_分析_Os02g基因鉴定/07_花青素位点SV/OsB2_Kala4/NIP_Kala4_locus.fasta
for tgt in "NIP:$NIPLOC" "MH63:MH63_Kala4_locus.fasta" "ZS97:ZS97_Kala4_locus.fasta"; do
  nm="${tgt%%:*}"; f="${tgt##*:}"
  makeblastdb -in "$f" -dbtype nucl -out tmpdb >/dev/null 2>&1
  echo "=== ZN65 OsB2 promoter-proximal vs $nm OsB2 locus ==="
  blastn -query osb2_prom_zn65.fa -db tmpdb -outfmt "6 qstart qend length pident" -evalue 1e-5 2>/dev/null | awk '{print "  qpos "$1"-"$2"  len="$3"  id="$4"%"}' | head -8
  cov=$(blastn -query osb2_prom_zn65.fa -db tmpdb -outfmt "6 length" -evalue 1e-5 2>/dev/null | awk '{s+=$1} END{print s+0}')
  echo "  -> total aligned bp of the 3.5kb promoter region: $cov"
done
echo DONE
