import sys
BP="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/V2/build_pkg.py"
s=open(BP,encoding="utf-8").read()
def rep(old,new,t):
    global s
    if s.count(old)!=1: print(f"ERR {t} count={s.count(old)}"); sys.exit(1)
    s=s.replace(old,new); print("ok",t)
# S1 caption
rep('P(d,"(a) Os02g0139500/OsOSC1 falls within the oxidosqualene-cyclase family (triterpene/sterol), not anthocyanin "\n  "enzymes. (b) The chromosome-2 MYB ZN652G3275 clusters with lignin/secondary-wall MYBs, not anthocyanin MYBs. "\n  "Maximum-likelihood trees (IQ-TREE).",size=9.5,color=GREY)',
'P(d,"(a) Os02g0139500/OsOSC1 falls within the 2,3-oxidosqualene-cyclase family (triterpene/sterol), not the "\n  "anthocyanin enzymes. (b) The chromosome-2 MYB ZN652G3275 clusters with lignin/secondary-wall MYBs (AtMYB58/63), "\n  "not anthocyanin MYBs. Methods: proteins aligned with MAFFT L-INS-i; maximum-likelihood trees inferred with "\n  "IQ-TREE 2 under a ModelFinder-selected model; node support is the ultrafast bootstrap (UFBoot, 1,000 replicates), "\n  "shown as percentages. Reference oxidosqualene cyclases: UniProt Q6Z2X6, P38605, Q2R712, H2KWF1, P48449. "\n  "Zero-length branches (e.g. ZN652G0345 / Os02g0139700) denote identical/collapsed sequences.",size=9.5,color=GREY)',"S1cap")
# S2 caption
rep('P(d,"Locus 98.93% identity, 64 SNPs; main structural difference is a 1.36-kb promoter TE insertion. CDS ~99.3% "\n  "identical. From the MAFFT locus alignment + exon coordinates.",size=9.5,color=GREY)',
'P(d,"OsOSC1 locus (Os02g0139500 / ZN652G0336): 98.93% identity over the aligned region, 64 SNPs, and a single "\n  "large structural difference \\u2014 a ~1.36-kb insertion in the ZN65 promoter annotated as a MULE-MuDR DNA "\n  "transposon (shown as a coloured block); CDS ~99.3% identical. Each indel is drawn to scale with its size and TE "\n  "family. Together with the OSC-family placement in Supplementary Fig. S1a, this confirms that the proposal gene "\n  "\\u201cOs02g\\u201d is an oxidosqualene cyclase, NOT an anthocyanin gene. From the MAFFT locus alignment and the "\n  "ZN65 exon coordinates.",size=9.5,color=GREY)',"S2cap")
open(BP,"w",encoding="utf-8").write(s)
print("SuppFig S1/S2 captions updated")
