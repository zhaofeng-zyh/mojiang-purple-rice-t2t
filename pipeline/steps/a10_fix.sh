source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
A10="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/A10"
GEN="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/genomes"
B="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定"
cd "$A10"
ZN65="/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa"
NIP="$B/08_全基因组SV_ZN65vsNIP/nip_genome.fa"
MH63="$B/10_多基因组Kala4/MH63RS3.fasta"; ZS97="$B/10_多基因组Kala4/ZS97RS3.fasta"
CEMPO="$GEN/CempoIreng_black.fna"; N22="$GEN/N22_aus.fna"; RUF="$GEN/Orufipogon_wild.fna"
PANEL="ZN65:$ZN65:japonica:purple Nipponbare:$NIP:japonica:white MH63:$MH63:indica:white ZS97:$ZS97:indica:white CempoIreng:$CEMPO:tropjaponica:black N22:$N22:aus:white Orufipogon:$RUF:wild:wild"
# conserved anchor = ZN65 locus 1-34131 (the part colinear with Nipponbare)
seqkit subseq -r "1:34131" zn65_locus.fa 2>/dev/null | seqkit replace -p '.*' -r 'anchor' 2>/dev/null > anchor.fa
echo -e "genome\tsubspecies\tpericarp\tortho_region(bp)\tinsert_localcov%\tinsert_localidy%\tcall" > PA_table_ortho.tsv
for e in $PANEL; do
  nm=${e%%:*}; r=${e#*:}; g=${r%%:*}; r=${r#*:}; ssp=${r%%:*}; per=${r##*:}
  # locate orthologous OsB2 locus via conserved anchor (single-copy)
  minimap2 -cx asm10 --secondary=no "$g" anchor.fa 2>/dev/null | sort -k11,11nr > an_$nm.paf
  read tname ts te < <(awk 'NR==1{print $6"\t"$8"\t"$9}' an_$nm.paf)
  if [ -z "$tname" ]; then echo -e "$nm\t$ssp\t$per\tNA\t0.0\t0.0\tNO_LOCUS" >> PA_table_ortho.tsv; continue; fi
  glen=$(awk -v c="$tname" '$1==c{print $2; exit}' "$g.fai" 2>/dev/null)
  [ -z "$glen" ] && { samtools faidx "$g" 2>/dev/null; glen=$(awk -v c="$tname" '$1==c{print $2; exit}' "$g.fai"); }
  lo=$(( ts-12000 > 1 ? ts-12000 : 1 )); hi=$(( te+12000 < glen ? te+12000 : glen ))
  samtools faidx "$g" "$tname:$lo-$hi" 2>/dev/null > ortho_$nm.fa
  makeblastdb -in ortho_$nm.fa -dbtype nucl -out odb_$nm >/dev/null 2>&1
  res=$(blastn -query insert.fa -db odb_$nm -evalue 1e-10 -outfmt '6 pident length qlen' 2>/dev/null | awk 'BEGIN{c=0;idy=0}{cov=$2/$3*100; if(cov>c){c=cov;idy=$1}}END{printf "%.1f\t%.1f",c,idy}')
  cov=$(echo "$res"|cut -f1)
  call=$(awk -v c="$cov" 'BEGIN{print (c>=60)?"PRESENT":((c>=20)?"PARTIAL":"ABSENT")}')
  reg=$((hi-lo))
  echo -e "$nm\t$ssp\t$per\t$reg\t$res\t$call" >> PA_table_ortho.tsv
  rm -f odb_$nm.*
done
echo "=== 正交(orthology-aware) 存在/缺失谱 ==="
column -t -s$'\t' PA_table_ortho.tsv
