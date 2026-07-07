import sys
M="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/manuscript.tex"
SI="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/supplementary/build_supplementary_tables.py"

def repl(s, old, new, tag):
    n=s.count(old)
    if n!=1:
        print(f"ERROR [{tag}]: anchor count = {n}"); sys.exit(1)
    return s.replace(old,new)

# ---------- manuscript.tex ----------
m=open(M,encoding="utf-8").read()

# (a) S8 nomenclature reference at first dual-name use
m=repl(m,
 "helix transcription factor (\\textit{OsB2}), respectively\\cite{oikawa2015,sun2018,kim2021}.",
 "helix transcription factor (\\textit{OsB2}), respectively\\cite{oikawa2015,sun2018,kim2021} "
 "(gene nomenclature for all pigmentation-related genes is summarised in Table~S8).",
 "S8ref")

# (b) pericarp vs hull distinction in Intro
m=repl(m,
 "of \\textit{OsB2} in the pericarp\\cite{oikawa2015,kim2021}.",
 "of \\textit{OsB2} in the pericarp\\cite{oikawa2015,kim2021}. Pericarp pigmentation "
 "(the \\textit{Kala1}/\\textit{Kala3}/\\textit{Kala4} system) is genetically distinct from "
 "hull and apiculus pigmentation, which is controlled by the \\textit{C}--\\textit{S}--\\textit{A} "
 "system\\cite{sun2018}; the present study concerns the pericarp.",
 "pericarp_hull")

# (c) OsOSC1 one-line framing in Discussion
m=repl(m,
 "into quantitatively analysable loci.",
 "into quantitatively analysable loci. We further note that a chromosome-2 locus "
 "(\\textit{Os02g0139500}) previously nominated as a pigmentation candidate is, by sequence "
 "alignment and phylogenetic analysis, an oxidosqualene cyclase (\\textit{OsOSC1}) of the "
 "triterpene/sterol pathway rather than an anthocyanin gene (Table~S7), underscoring the "
 "importance of resolving candidate-gene identity before functional interpretation.",
 "OsOSC1line")

open(M,"w",encoding="utf-8").write(m)
print("manuscript.tex: 3 edits OK")

# ---------- build_supplementary_tables.py ----------
b=open(SI,encoding="utf-8").read()

# (d) fix S6 dating inconsistency 1.4 -> 0.2 Mya
b=repl(b,
 "complete element; LTRs dated ~1.4 Mya",
 "complete element; LTRs dated ~0.2 Mya (K2P; 95% CI 0.04-0.39 Myr)",
 "S6date")

# (e) add S8 (nomenclature) + S9 (interim polyphenol) before output line
s89='''# --- S8 Dual nomenclature ---
s8 = [["OsC1", "Kala3", "Os06g0205100", "ZN656G0716", "R2R3-MYB activator (MBW)"],
      ["OsB2", "Kala4", "Os04g0557500", "ZN654G2687", "bHLH activator (MBW); anthocyanin"],
      ["OsDFR", "Kala1", "Os01g0633500", "ZN651G2772", "dihydroflavonol-4-reductase (structural)"],
      ["Rc", "-", "Os07g0211500", "ZN657G0823", "bHLH; proanthocyanidin regulator"],
      ["OsTTG1", "-", "Os02g0682500", "ZN652G3195", "WD40 partner (MBW)"],
      ["OsOSC1", "proposal 'Os02g'", "Os02g0139500", "ZN652G0336", "2,3-oxidosqualene cyclase (triterpene/sterol; NOT anthocyanin)"]]
sheet("S8_Nomenclature", "Table S8. Gene nomenclature used in this study: common/literature symbol, Kala designation, RAP-DB (Nipponbare) locus and ZN65 gene ID",
      ["Symbol", "Kala designation", "RAP-DB (Nipponbare)", "ZN65 gene", "Role"], s8)

# --- S9 Interim total-polyphenol (NOT anthocyanin-specific) ---
s9 = [["Black rice", "39.65 +/- 0.18", "Folin-type kit, gallic-acid equiv., ethanol-water extract"],
      ["Purple rice (Mojiang, ZN65 type)", "25.02 +/- 0.15", "same"],
      ["Red rice", "24.65 +/- 0.19", "same"]]
sheet("S9_TotalPolyphenol_interim", "Table S9. Prior in-house total-polyphenol comparison (interim context; Zhang 2024 thesis). NOTE: total polyphenol on an extract basis (ug/mL), NOT anthocyanin-specific and NOT a grain/pericarp content; HPLC-DAD pericarp anthocyanin quantification is in progress.",
      ["Sample", "Total polyphenol (ug/mL extract)", "Method"], s9)

'''
b=repl(b,
 'out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Supplementary_Tables.xlsx")',
 s89+'out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Supplementary_Tables.xlsx")',
 "S8S9insert")

open(SI,"w",encoding="utf-8").write(b)
print("build_supplementary_tables.py: S6 fix + S8/S9 added OK")
