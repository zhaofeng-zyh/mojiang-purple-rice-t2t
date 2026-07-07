f="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/supplementary/build_supplementary_tables.py"
s=open(f).read()
old='["LTR Assembly Index (LAI)", "23.3", "reference-grade (>20)"],'
new='["LTR Assembly Index (LAI)", "15.5", "reference-quality range 10-20 (Ou 2018); LTR_retriever"],'
if old in s: open(f,"w").write(s.replace(old,new)); print("S1 build script LAI fixed")
else: print("S1 old string not found")
