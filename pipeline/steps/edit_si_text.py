import sys
SI="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/supplementary/supplementary_information.tex"
def repl(s,old,new,tag):
    n=s.count(old)
    if n!=1: print(f"ERROR [{tag}] count={n}"); sys.exit(1)
    return s.replace(old,new)
s=open(SI,encoding="utf-8").read()

# (1) FIX dating contradiction: K=0.036 (~1.4 Mya) -> K2P ~0.2 Mya + CI
old1=("retained both long terminal repeats (5$'$ 994~bp; 3$'$ 1{,}026~bp). The LTRs were aligned\n"
"globally (EMBOSS \\texttt{stretcher}); their divergence $K=0.036$ was converted to an age\n"
"$T=K/(2\\mu)$ with $\\mu=1.3\\times10^{-8}$ substitutions per site per year.")
new1=("retained both long terminal repeats. The two LTRs were aligned (MAFFT; raw indels\n"
"excluded, 994 aligned sites) and differed by only $\\sim$5 substitutions; the\n"
"Kimura-2-parameter distance $K$ was converted to an age $T=K/(2\\mu)$ with\n"
"$\\mu=1.3\\times10^{-8}$ substitutions per site per year, giving $\\approx$0.2 million\n"
"years (95\\% bootstrap CI 0.04--0.39 Myr; 1{,}000 column-bootstrap replicates). This\n"
"young intronic element predates rice domestication and is therefore unlikely to be the\n"
"pigmentation-causing change itself.")
s=repl(s,old1,new1,"dating")

# (2) expand comparative-genomics panel + add A10 orthology test + A12 3K projection
old2=("Reference genomes used\n"
"were Nipponbare (IRGSP-1.0; Ensembl Plants release-58), the gap-free indica references\n"
"MH63RS3 and ZS97RS3 (RIGW/HZAU), and the black-pericarp landrace Cempo Ireng\n"
"(GenBank GCA\\_055776245.1).")
new2=("Reference genomes used\n"
"were Nipponbare (IRGSP-1.0; Ensembl Plants release-58), the gap-free indica references\n"
"MH63RS3 and ZS97RS3 (RIGW/HZAU), the black-pericarp landrace Cempo Ireng\n"
"(GenBank GCA\\_055776245.1), the aus accession N22 (GCA\\_001952365.2) and the wild\n"
"progenitor \\textit{Oryza rufipogon} W1943 (GCA\\_000817225.1).")
s=repl(s,old2,new2,"panel")

# (3) add two methods paragraphs (OsB2 insertion P/A; WD40; subspecies 3K) before Expression
anchor="\\textbf{Expression.}"
addnew=("\\textbf{OsB2-promoter insertion presence/absence and WD40 partner.} The 5.9-kb\n"
"ZN65-specific \\textit{OsB2} distal-promoter segment (the portion of the ZN65 locus with\n"
"no aligned Nipponbare counterpart) was tested for presence at the orthologous locus in\n"
"each of the seven panel genomes by anchoring on the 34-kb conserved \\textit{OsB2} region\n"
"and scoring the insertion only within that orthologous interval (Table~S10). The MBW\n"
"WD40 partner was identified by BLAST of Arabidopsis TTG1 (UniProt Q9XGN1) against the\n"
"ZN65 proteome and confirmed by reciprocal best hits against Nipponbare; ZN65\n"
"OsTTG1=ZN652G3195 (Table~S5).\n\n"
"\\textbf{Subspecies placement.} ZN65 was genotyped at the 1{,}011{,}601 pruned SNPs of the\n"
"3{,}000 Rice Genomes Project (v2.1) from its whole-genome alignment to Nipponbare, merged\n"
"into the 3{,}024-accession panel (PLINK~1.9) and projected by principal-component analysis\n"
"(PLINK~2.0); all 50 nearest accessions, and the nearest subpopulation centroid, were\n"
"\\textit{indica} (Table~S11).\n\n")
s=repl(s,anchor,addnew+anchor,"newmethods")
open(SI,"w",encoding="utf-8").write(s)
print("supplementary_information.tex: dating fixed + panel/A10/A12 added")
