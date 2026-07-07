f="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/manuscript.tex"
s=open(f).read()
edits=[
# R1 sentence
("""an LTR Assembly Index of 23.3 (above the
reference-quality threshold of 20), read-mapping rates of 99.96\\% (HiFi), 99.99\\%""",
"""an LTR Assembly Index of 15.5 (LTR\\_retriever;
within the reference-quality range of 10--20)\\cite{ou2018lai}, read-mapping rates of 99.96\\% (HiFi), 99.99\\%"""),
# Methods mention
("Assembly Index (LAI~23.3), read-mapping rates",
 "Assembly Index (LAI~15.5, computed with LTR\\_retriever\\cite{ou2018lai}), read-mapping rates"),
# add bibitem before end of bibliography
("\\end{thebibliography}",
 """\\bibitem{ou2018lai} Ou, S., Chen, J. \\& Jiang, N. Assessing genome assembly quality using the LTR Assembly Index (LAI). \\textit{Nucleic Acids Res.} \\textbf{46}, e126 (2018). doi:10.1093/nar/gky730.
\\end{thebibliography}"""),
]
n=0
for old,new in edits:
    if old in s: s=s.replace(old,new,1); n+=1
    else: print("NOT FOUND:", old[:50])
open(f,"w").write(s)
print(f"applied {n}/3 edits")
