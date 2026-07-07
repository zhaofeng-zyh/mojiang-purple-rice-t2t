#!/bin/bash
# 从当前 .tex 重新生成 Word 版 (图用 png 以便 Word 显示)。Mac 上跑 (pandoc 在 zn65phylo env)。
PANDOC="/opt/homebrew/Caskroom/miniconda/base/envs/zn65phylo/bin/pandoc"
D="/Volumes/data2/墨江紫米研究/12_论文Paper1_Manuscript"; FIG="$D/figures"
cd "$D"
# 主文 -> docx (pdf图改png)
python3 -c "
s=open('manuscript.tex',encoding='utf-8').read()
for fn in ['Fig1_genome_landscape','Fig_pathway_complement','Fig3_pigmentation_loci','Fig4_kala4_architecture']: s=s.replace(fn+'.pdf',fn+'.png')
open('manuscript_docx.tex','w',encoding='utf-8').write(s)"
"$PANDOC" manuscript_docx.tex -o manuscript.docx --resource-path="$D:$FIG"
# 补充材料 -> docx
cd "$D/supplementary"
python3 -c "
s=open('supplementary_information.tex',encoding='utf-8').read()
for fn in ['SuppFig_Os02g_phylogenies','SuppFig_OsOSC1_variants']: s=s.replace(fn+'.pdf',fn+'.png')
open('supplementary_information_docx.tex','w',encoding='utf-8').write(s)"
"$PANDOC" supplementary_information_docx.tex -o Supplementary_Information.docx --resource-path="$D/supplementary:$FIG"
echo "Word 版已重建: manuscript.docx + Supplementary_Information.docx"
