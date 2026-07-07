f="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/manuscript.tex"
s=open(f).read()
old="""every chromosome, with a cumulative excess of roughly 22\\,Mb attributable to
transposon expansion. The most prominent rearrangement was a large inversion on"""
new="""every chromosome, with a cumulative excess of roughly 22\\,Mb attributable to
transposon expansion. Of the 73\\,Mb of ZN65 sequence with no aligned Nipponbare
counterpart, 45\\,Mb (62\\%) was repetitive and dominated by LTR/\\textit{Gypsy}
retrotransposons (26.9\\,Mb), with DNA transposons contributing a further 11.1\\,Mb,
identifying \\textit{Gypsy} proliferation as the principal driver of the
lineage-specific sequence gain. The most prominent rearrangement was a large inversion on"""
if old in s: open(f,"w").write(s.replace(old,new,1)); print("A16 R2 edit OK")
else: print("NOT FOUND")
