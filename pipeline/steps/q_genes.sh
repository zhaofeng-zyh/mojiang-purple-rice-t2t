#!/bin/bash
M=/mnt/data2/墨江紫米研究/14_转录组_重比对ZN65T2T/counts/ZN65_gene_count_matrix.tsv
echo -e "ZN65_gene\tlocus\tSD1\tSD2\tSD3\tTG1\tTG3"
declare -A name=( [ZN654G2687]="OsB2/Kala4" [ZN654G2685]="OsB1" [ZN656G0716]="OsC1/Kala3" [ZN651G2772]="Kala1/OsDFR" [ZN657G0823]="Rc" [ZN652G0336]="OsOSC1" )
for g in ZN654G2687 ZN654G2685 ZN656G0716 ZN651G2772 ZN657G0823 ZN652G0336; do
  row=$(grep -P "^${g}\t" "$M")
  echo -e "${g}\t${name[$g]}\t$(echo "$row" | cut -f2-6)"
done
