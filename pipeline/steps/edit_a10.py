import sys
M="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/manuscript.tex"
SI="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/supplementary/build_supplementary_tables.py"
def repl(s,old,new,tag):
    n=s.count(old)
    if n!=1: print(f"ERROR [{tag}] count={n}"); sys.exit(1)
    return s.replace(old,new)
m=open(M,encoding="utf-8").read()

# Edit 1: Results — 7-genome orthology-anchored P/A
old1=("place this configuration among other rice accessions, we\n"
"compared the \\textit{Kala4}/\\textit{OsB2} locus across ZN65, the white-pericarp\n"
"references Nipponbare, MH63 and ZS97, and the black-pericarp landrace Cempo Ireng.\n"
"The proximal $\\sim$1\\,kb of the OsB2 promoter was conserved among the japonica-type\n"
"genomes ZN65, Nipponbare and Cempo Ireng ($>$98\\% identity), but the distal\n"
"transposon-containing promoter segment did not align between ZN65 and any other\n"
"genome, and the locus as a whole was structurally hypervariable, aligning over only\n"
"39\\% of its length between ZN65 and Cempo Ireng despite 97\\% identity in the aligned\n"
"portion.")
new1=("test whether this configuration is shared or lineage-specific, we examined the\n"
"\\textit{Kala4}/\\textit{OsB2} locus across a panel of seven genomes spanning the major\n"
"rice varietal groups and pericarp colours: ZN65 (purple), the white-pericarp references\n"
"Nipponbare (temperate japonica), MH63 and ZS97 (indica) and N22 (aus), the\n"
"black-pericarp landrace Cempo Ireng (tropical japonica), and the wild progenitor\n"
"\\textit{Oryza rufipogon} (W1943). After anchoring on the 34-kb conserved \\textit{OsB2}\n"
"region to locate the orthologous locus in each genome, we tested whether the\n"
"ZN65-specific 5.9-kb distal-promoter insertion was present there. The insertion was\n"
"detected only in ZN65 and was absent from the orthologous \\textit{OsB2} locus in all\n"
"six other genomes---including the black-pericarp Cempo Ireng and the wild\n"
"progenitor---whereas the proximal $\\sim$1\\,kb of the promoter remained conserved across\n"
"the panel ($>$98\\% identity) (Table~S10). The distal transposon-containing segment\n"
"therefore represents a ZN65-specific insertion rather than a feature shared with\n"
"canonical pigmented or ancestral alleles.")
m=repl(m,old1,new1,"Results_A10")

# Edit 2: Discussion — panel support replaces n=1
old2=("as shown by the extensive structural divergence between ZN65 and the\n"
"black landrace Cempo Ireng.")
new2=("as shown by the absence of the ZN65-specific OsB2-promoter insertion at the\n"
"orthologous locus across a seven-genome panel (japonica, indica, aus and the wild\n"
"progenitor \\textit{O. rufipogon}), including the black landrace Cempo Ireng (Table~S10).")
m=repl(m,old2,new2,"Discussion_A10")
open(M,"w",encoding="utf-8").write(m)
print("manuscript.tex: A10 Results+Discussion OK")

# Edit 3: SI Table S10
b=open(SI,encoding="utf-8").read()
s10code='''# --- S10 multi-genome OsB2 insertion presence/absence (A10) ---
s10=[["ZN65","purple-pericarp landrace (this study)","PRESENT (100% local coverage)"],
     ["Nipponbare","temperate japonica / white","absent"],
     ["MH63","indica / white","absent"],
     ["ZS97","indica / white","absent"],
     ["Cempo Ireng","tropical japonica / black","absent"],
     ["N22","aus / white","absent"],
     ["O. rufipogon (W1943)","wild progenitor","absent"]]
sheet("S10_OsB2_insertion_PA","Table S10. Presence/absence of the ZN65-specific 5.9-kb OsB2 distal-promoter insertion at the ORTHOLOGOUS Kala4/OsB2 locus across a seven-genome rice panel (orthology-anchored test; insertion is ZN65-lineage-specific, absent even from black rice and the wild progenitor)",["Genome","Group / pericarp","ZN65-specific insertion at orthologous locus"],s10)

'''
b=repl(b,'out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Supplementary_Tables.xlsx")',
       s10code+'out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Supplementary_Tables.xlsx")',"S10")
open(SI,"w",encoding="utf-8").write(b)
print("build_supplementary_tables.py: S10 added OK")
