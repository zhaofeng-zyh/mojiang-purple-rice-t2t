import sys
M="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/manuscript.tex"
s=open(M,encoding="utf-8").read()
old=("Orthologues were identified by reciprocal BLAST against\n"
"the IRGSP-1.0 proteome.")
new=("Orthologues were identified by reciprocal BLAST against\n"
"the IRGSP-1.0 proteome. Presence or absence of the ZN65-specific \\textit{OsB2}\n"
"distal-promoter insertion was assessed across a seven-genome panel that additionally\n"
"included N22 (aus; GenBank GCA\\_001952365.2) and \\textit{Oryza rufipogon} W1943\n"
"(wild progenitor; GCA\\_000817225.1); the orthologous \\textit{OsB2} locus was located\n"
"in each genome by anchoring on the 34-kb conserved region, and the 5.9-kb insertion\n"
"was scored only within the orthologous interval to avoid false positives from\n"
"dispersed transposon copies.")
n=s.count(old)
if n!=1: print("ERROR count",n); sys.exit(1)
open(M,"w",encoding="utf-8").write(s.replace(old,new))
print("Methods: A10 panel + orthology-anchored test added OK")
