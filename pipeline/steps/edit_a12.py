import sys
M="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/manuscript.tex"
SI="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/supplementary/build_supplementary_tables.py"
def repl(s,old,new,tag):
    n=s.count(old)
    if n!=1: print(f"ERROR [{tag}] count={n}"); sys.exit(1)
    return s.replace(old,new)
m=open(M,encoding="utf-8").read()

# (1) R2: japonica -> indica via 3K projection
old1=("The substitution density and locus-level comparisons (below) placed\n"
"ZN65 within the japonica clade rather than with the indica accessions Minghui 63 and\n"
"Zhenshan 97.")
new1=("To place ZN65 within \\textit{Oryza sativa}, we projected its genotypes at 1.0\\,million\n"
"pruned SNPs onto the 3,000 Rice Genomes Project panel (3,024 accessions with nine\n"
"subpopulation assignments). ZN65 fell within the \\textit{indica} group: all 50 nearest\n"
"accessions were \\textit{indica}, and ZN65 lay closest to the \\textit{indica}\n"
"subpopulation centroid, 2.5-fold nearer than to any \\textit{japonica}, \\textit{aus} or\n"
"aromatic group (Table~S11). ZN65 is therefore an \\textit{indica} landrace; notably, its\n"
"pigmentation-locus haplotypes (below) are more conserved with the japonica reference\n"
"than its genome-wide indica background would predict.")
m=repl(m,old1,new1,"R2_indica")

# (2) Methods: add 3K projection
old2="Gene expression was quantified from seedling and tillering RNA-seq."
new2=("Gene expression was quantified from seedling and tillering RNA-seq. Subspecies\n"
"placement used the 3,000 Rice Genomes Project pruned SNP set (v2.1;\n"
"1,011,601 SNPs)\\cite{wang20183k}: ZN65 was genotyped at these positions from its\n"
"whole-genome alignment to Nipponbare, merged into the 3,024-accession panel with\n"
"PLINK\\,1.9, and projected by principal-component analysis (PLINK\\,2.0); subpopulation\n"
"labels were taken from the published K=9 admixture matrix.")
m=repl(m,old2,new2,"Methods_3K")
open(M,"w",encoding="utf-8").write(m)
print("manuscript.tex: A12 indica correction + Methods OK")

# (3) bib entry for 3K project
b=open(M,encoding="utf-8").read()
if "wang20183k" not in b.split("\\bibitem{wang20183k}")[0] or b.count("\\bibitem{wang20183k}")==0:
    anchor="\\bibitem{nattestad2016assemblytics}"
    add=("\\bibitem{wang20183k} Wang, W. \\textit{et al.} Genomic variation in 3,010 diverse\n"
"accessions of Asian cultivated rice. \\textit{Nature} \\textbf{557}, 43--49 (2018).\n"
"doi:10.1038/s41586-018-0063-9.\n")
    b=repl(b,anchor,add+anchor,"bib3k")
    open(M,"w",encoding="utf-8").write(b)
    print("bib: 3K Project citation added")

# (4) SI Table S11
si=open(SI,encoding="utf-8").read()
s11='''# --- S11 ZN65 subspecies assignment (3K-RGP projection, A12) ---
s11=[["Nearest 50 accessions","50/50 indica","K=9 subpopulation labels (Qmatrix-k9-3kRG)"],
     ["Distance to indica centroid","0.0194","closest"],
     ["Distance to tropical-japonica centroid","0.0497","2.6x farther"],
     ["Distance to temperate-japonica centroid","0.0581","3.0x farther"],
     ["Distance to aus centroid","0.0686","3.5x farther"],
     ["Distance to aromatic centroid","0.1014","5.2x farther"],
     ["Assignment","indica","projected onto 3,024 3K-RGP accessions, 1.01M pruned SNPs, PLINK2 PCA"]]
sheet("S11_ZN65_subspecies","Table S11. Subspecies placement of ZN65 by projection onto the 3,000 Rice Genomes Project panel (ZN65 = indica)",["Metric","Value","Note"],s11)

'''
si=repl(si,'out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Supplementary_Tables.xlsx")',
        s11+'out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Supplementary_Tables.xlsx")',"S11")
open(SI,"w",encoding="utf-8").write(si)
print("SI: Table S11 added")
