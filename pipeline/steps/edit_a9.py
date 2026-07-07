f="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/manuscript.tex"
s=open(f).read()
e=[
("summary correspondingly groups these, e.g.\\ 449 translocations).",
 """summary correspondingly groups these, e.g.\\ 449 translocations). Robustness of the
structural-variant calls to alignment stringency was confirmed by repeating the
ZN65--Nipponbare comparison at the asm10 and asm20 minimap2 presets, which recovered
comparable structural counts (121, 132 and 126 inversions, and 432, 427 and 399
duplications, at asm5, asm10 and asm20, respectively), and by an orthogonal
nucmer/Assemblytics analysis\\cite{nattestad2016assemblytics}."""),
("\\end{thebibliography}",
 """\\bibitem{nattestad2016assemblytics} Nattestad, M. \\& Schatz, M.C. Assemblytics: a web analytics tool for the detection of variants from an assembly. \\textit{Bioinformatics} \\textbf{32}, 3021--3023 (2016). doi:10.1093/bioinformatics/btw369.
\\end{thebibliography}"""),
]
n=0
for o,nw in e:
    if o in s: s=s.replace(o,nw,1); n+=1
    else: print("NOTFOUND:",o[:40])
open(f,"w").write(s); print(f"{n}/2 applied")
