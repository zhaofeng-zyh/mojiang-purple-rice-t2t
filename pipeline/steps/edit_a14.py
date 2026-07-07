import sys
M="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/manuscript.tex"
s=open(M,encoding="utf-8").read()
old="otherwise being 98.9\\% identical. ZN65\ntherefore carries the complete functional regulatory and structural allele complement"
ins=("otherwise being 98.9\\% identical. Completing the MYB--bHLH--WD40 module, ZN65 also\n"
     "encodes a full-length WD40 partner: its \\textit{OsTTG1} orthologue (ZN652G3195,\n"
     "355 residues, intact open reading frame) is 62\\% identical to Arabidopsis TTG1 and\n"
     "99.7\\% identical to the Nipponbare \\textit{OsTTG1} (\\textit{Os02g0682500}), and is\n"
     "constitutively expressed across vegetative tissues, so that all three MBW\n"
     "components---the R2R3-MYB \\textit{OsC1}, the bHLH \\textit{OsB2} and the WD40\n"
     "\\textit{OsTTG1}---are present and intact in ZN65. ZN65\n"
     "therefore carries the complete functional regulatory and structural allele complement")
n=s.count(old)
if n!=1:
    print("ERROR: anchor count =",n); sys.exit(1)
open(M,"w",encoding="utf-8").write(s.replace(old,ins))
print("OK: inserted A14 WD40/OsTTG1 sentence (1 substitution)")
