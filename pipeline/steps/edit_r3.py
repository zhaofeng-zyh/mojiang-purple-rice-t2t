f="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/manuscript.tex"
s=open(f).read()
old="""providing a
complete enzymatic complement for pericarp anthocyanin biosynthesis; copy numbers
relative to Nipponbare were not formally tested here."""
new="""providing a
complete enzymatic complement for pericarp anthocyanin biosynthesis. All 100 ZN65
flavonoid-pathway genes had high-identity orthologues in Nipponbare (mean amino-acid
identity 84--100\\% per enzyme step), indicating that the pathway is structurally
conserved between the two genomes rather than expanded or eroded in ZN65; the
pigmentation difference therefore maps to the regulatory loci rather than to the
biosynthetic enzymes themselves (Table~S4)."""
if old in s:
    open(f,"w").write(s.replace(old,new)); print("R3 edited OK")
else:
    print("OLD STRING NOT FOUND")
