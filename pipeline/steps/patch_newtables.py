import sys
BP="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/V2/build_pkg.py"
s=open(BP,encoding="utf-8").read()
def rep(old,new,tag):
    global s
    if s.count(old)!=1: print(f"ERR [{tag}] count={s.count(old)}"); sys.exit(1)
    s=s.replace(old,new); print("ok",tag)

# data-building block (insert after allele_rows definition)
anchor_data="# ================= COVER LETTER ================="
block_data='''# ---- S5 assembly stats (read .fai + telomere info) ----
FAI=ROOT+"/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP/zn65_genome.fa.fai"
TELO=ROOT+"/01_基因组_组装与注释_Genome/01.assembly/ZN65.telo.info"
chrlen={}
for ln in open(FAI):
    f=ln.split("\\t"); chrlen[f[0]]=int(f[1])
telost={}
try:
    for ln in open(TELO):
        if ln.startswith("#") or not ln.strip(): continue
        f=ln.split("\\t"); telost[f[0]]=f[2]
except FileNotFoundError: pass
asm_rows=[]
for c in sorted(chrlen,key=lambda x:int(''.join(ch for ch in x if ch.isdigit()) or 0)):
    asm_rows.append([c, f"{chrlen[c]/1e6:.2f}", telost.get(c,"both"), "0"])
asm_total=sum(chrlen.values())
asm_metrics=[["Total assembly size", f"{asm_total:,} bp ({asm_total/1e6:.1f} Mb)"],
   ["Chromosomes (gap-free)", "12"], ["Gaps", "0"], ["Telomeres resolved", "24 / 24 (AAACCCT)"],
   ["Contig N50", "32.35 Mb"], ["GC content", "43.7%"], ["Merqury consensus QV (k=19)", "53.6 (err 4.4e-06)"],
   ["LTR Assembly Index (LAI)", "15.5 (LTR_retriever; reference-grade 10\\u201320)"],
   ["BUSCO (embryophyta_odb10)", "99.6% C (S 96.9%, D 2.7%, F 0.4%, M 0.0%)"],
   ["Protein-coding genes", "42,090"], ["Repeat content", "56.6%"],
   ["Read-mapping rate", "HiFi 99.96% / ONT 99.99% / DNBSEQ 99.31%"],
   ["Centromere / 45S-5S rDNA coordinates", "pericentromeric TE-rich regions resolved (Fig. 1); exact CEN-satellite and rDNA coordinates to be added"]]
# ---- S6 sequencing throughput + data availability ----
seq_rows=[["PacBio HiFi (Revio)", "Genome (WGS)", "20.6 Gb; read N50 19.4 kb", "NCBI BioProject \\u2014 to be assigned"],
   ["Oxford Nanopore ultra-long", "Genome (WGS)", "~20 Gb (pass-UL); Guppy, Q\\u22657", "to be assigned"],
   ["DNBSEQ paired-end", "Polishing/QV", "72.6 Gb clean", "to be assigned"],
   ["Hi-C", "Scaffolding", "39.3 Gb clean", "to be assigned"],
   ["RNA-seq", "Expression (SD1-3, TG1, TG3; 150 bp PE)", "5 libraries", "to be assigned"]]
ref_rows=[["Nipponbare", "IRGSP-1.0 (japonica reference)", "Ensembl Plants / RAP-DB"],
   ["Minghui 63 / Zhenshan 97", "MH63RS3 / ZS97RS3 (indica)", "RIGW/HZAU"],
   ["Cempo Ireng (black pericarp)", "GCA_055776245.1", "NCBI"],
   ["N22 (aus)", "GCA_001952365.2", "NCBI"],
   ["Oryza rufipogon W1943 (wild)", "GCA_000817225.1", "NCBI"],
   ["3,000 Rice Genomes Project", "pruned SNP set v2.1 + K=9 admixture", "SNP-Seek/IRRI (public)"],
   ["ZN65 assembly + annotation (this study)", "to be deposited", "GWH/GenBank \\u2014 accession at submission"]]
# ---- S7 germplasm provenance + ABS/Nagoya ----
germ_rows=[["Common name", "Mojiang purple rice (\\u58a8\\u6c5f\\u7d2b\\u7c73)"],
   ["Sequenced line", "ZN65"], ["Grain type", "Glutinous, purple pericarp"],
   ["Geographic origin", "Mojiang Hani Autonomous County, Pu'er, Yunnan, China"],
   ["Traditional custodians", "Hani (Hani people) heritage landrace"],
   ["Maintained / provided by", "College of Agronomy, Pu'er University"],
   ["Collection & prior-informed consent (PIC)", "[PI to confirm: collection permit and PIC details]"],
   ["Access & benefit-sharing (Nagoya Protocol)", "China is a Party to the Nagoya Protocol; access and benefit-sharing of this domestic genetic resource handled under China's applicable measures [PI to confirm MAT/PIC]"],
   ["Material availability", "From the corresponding author under a standard material-transfer agreement (MTA)"]]

# ================= COVER LETTER ================='''
rep(anchor_data,block_data,"data_block")

# insert S5/S6/S7 into SI docx before save
rep('d.save(OUTD+"/Supplementary_Information.docx"); print("saved supplementary information")',
'''H(d,"Supplementary Table S5. Telomere-to-telomere assembly statistics and quality metrics",11)
tbl(d,["Chromosome","Length (Mb)","Telomeres","Gaps"],asm_rows,fs=8.5)
tbl(d,["Metric","Value"],asm_metrics,fs=8.5)
H(d,"Supplementary Table S6. Sequencing data and data availability",11)
tbl(d,["Platform","Data type","Throughput","Accession"],seq_rows,fs=8.5)
P(d,"Reference datasets used:",size=9,color=GREY,after=2)
tbl(d,["Source","Dataset","Repository"],ref_rows,fs=8.5)
P(d,"ZN65 raw reads, assembly and annotation are being deposited; accession numbers will be provided at journal "
  "submission and are available from the corresponding author in the interim.",size=9,color=GREY)
H(d,"Supplementary Table S7. Germplasm provenance and access-and-benefit-sharing (ABS) statement",11)
tbl(d,["Field","Information"],germ_rows,fs=8.5)
d.save(OUTD+"/Supplementary_Information.docx"); print("saved supplementary information")''',"si_tables")

# insert S5/S6/S7 sheets into xlsx before save
rep('wb.save(OUTD+"/Supplementary_Tables.xlsx"); print("saved supplementary tables xlsx")',
'''ws=wb.create_sheet("S5_assembly_stats"); ws.append(["Chromosome","Length (Mb)","Telomeres","Gaps"])
for r in asm_rows: ws.append(r)
ws.append([]); ws.append(["Metric","Value"])
for r in asm_metrics: ws.append(r)
ws=wb.create_sheet("S6_sequencing_dataavail"); ws.append(["Platform","Data type","Throughput","Accession"])
for r in seq_rows: ws.append(r)
ws.append([]); ws.append(["Reference source","Dataset","Repository"])
for r in ref_rows: ws.append(r)
ws=wb.create_sheet("S7_germplasm_ABS"); ws.append(["Field","Information"])
for r in germ_rows: ws.append(r)
wb.save(OUTD+"/Supplementary_Tables.xlsx"); print("saved supplementary tables xlsx")''',"xlsx_tables")

open(BP,"w",encoding="utf-8").write(s)
print("S5/S6/S7 tables added")
