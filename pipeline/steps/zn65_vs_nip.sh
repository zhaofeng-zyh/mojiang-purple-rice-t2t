source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
WD="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/23_亚种归属_3KRG_A12"
mkdir -p "$WD/zn65gt"
ZN65="/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa"
NIP="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP/nip_genome.fa"
echo "ZN65->NIP asm10 start $(date)"
minimap2 -cx asm10 --cs -t 20 "$NIP" "$ZN65" 2>/dev/null | sort -k6,6 -k8,8n 2>/dev/null | paftools.js call -L5000 - 2>/dev/null > "$WD/zn65gt/zn65_vs_nip.var"
echo "done $(date); R+V lines:"; cut -f1 "$WD/zn65gt/zn65_vs_nip.var" | sort | uniq -c
touch "$WD/zn65gt/ALIGN_DONE"
