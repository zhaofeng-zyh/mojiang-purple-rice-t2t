import sys
SI="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/supplementary/build_supplementary_tables.py"
s=open(SI,encoding="utf-8").read()
anchor='"~19% non-aligned; promoter Stowaway/MuDR"]]'
row=('"~19% non-aligned; promoter Stowaway/MuDR"],\n'
     '      ["OsTTG1 (MBW WD40 partner)", "Os02g0682500", "ZN652G3195", "Chr2", "99.72", "355", "%NIPAA%", '
     '"Functional (full-length 355 aa, intact ORF; 62% id to Arabidopsis TTG1; reciprocal best hit both ways)", '
     '"Conserved (99.7% id to Nipponbare); constitutively expressed in vegetative tissue (Step-14 RNA-seq)"]]')
n=s.count(anchor)
if n!=1:
    print("ERROR anchor count",n); sys.exit(1)
s=s.replace(anchor,row)
open(SI,"w",encoding="utf-8").write(s)
print("OK: S5 OsTTG1 row added (placeholder %NIPAA% for Nipponbare aa)")
