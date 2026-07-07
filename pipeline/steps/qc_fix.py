import sys
M="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/manuscript.tex"
def repl(s,old,new,tag):
    n=s.count(old)
    if n!=1: print(f"ERROR [{tag}] count={n}"); sys.exit(1)
    return s.replace(old,new)
s=open(M,encoding="utf-8").read()

# 1) FIX title contradiction
s=repl(s,
 "\\subsection{The flavonoid pathway is complete and copy-number-expanded}",
 "\\subsection{The flavonoid pathway is complete and structurally conserved}","title")

# 2) cite Table S2 (annotation)
s=repl(s,
 "Annotation identified 42,090 protein-coding\ngenes, 1,099 microRNAs and additional non-coding RNAs, while repetitive elements",
 "Annotation identified 42,090 protein-coding\ngenes, 1,099 microRNAs and additional non-coding RNAs (Table~S2), while repetitive elements","S2")

# 3) cite Table S3 (genome-wide SV)
s=repl(s,
 "432 duplications and\n357 inverted duplications (Fig.~\\ref{fig:sv}).",
 "432 duplications and\n357 inverted duplications (Fig.~\\ref{fig:sv}; Table~S3).","S3")

# 4) cite Table S6 (OsB2 TEs)
s=repl(s,
 "the proximal promoter carried a nested\ncluster of SINE, Helitron, hAT and PIF-Harbinger elements.",
 "the proximal promoter carried a nested\ncluster of SINE, Helitron, hAT and PIF-Harbinger elements (Table~S6).","S6")

# 5) cite Table S9 (interim polyphenol) in Discussion
s=repl(s,
 "pericarp developmental transcriptomes and targeted\nanthocyanin quantification are the natural next step.",
 "earlier in-house assays found higher total\npolyphenol in pigmented than in white grain (Table~S9), but locus-resolved anthocyanin\nquantification by HPLC-DAD remains the natural next step alongside pericarp\ndevelopmental transcriptomes.","S9")

# 6) abstract: state indica placement
s=repl(s,
 "every chromosome, reflecting transposon expansion. Systematic dissection of the",
 "every chromosome, reflecting transposon expansion. Projection onto the 3,000 Rice\nGenomes panel placed ZN65 within the \\textit{indica} group. Systematic dissection of the","abs_indica")

# 7) improve Data availability (concrete, preprint-grade)
s=repl(s,
 "Genome assembly, annotation and raw reads will be deposited in NCBI/GSA under a\nBioProject accession upon submission.",
 "The raw sequencing data (PacBio HiFi, Oxford Nanopore ultra-long, Hi-C, DNBSEQ and\nRNA-seq), the telomere-to-telomere genome assembly and its annotation are being\ndeposited in the NCBI BioProject / Genome Sequence Archive (GSA); accession numbers\nwill be added before journal submission and are available from the authors on request in\nthe interim. The 3,000 Rice Genomes Project genotype data used for subspecies placement\nare publicly available\\cite{wang20183k}; all other reference genomes are available under\nthe accessions cited in Methods.","dataavail")

open(M,"w",encoding="utf-8").write(s)
print("QC fixes applied: title, S2/S3/S6/S9 cites, abstract indica, data-availability")
