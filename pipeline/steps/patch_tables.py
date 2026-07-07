import sys
BP="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/V2/build_pkg.py"
s=open(BP,encoding="utf-8").read()
def rep(old,new,tag):
    global s
    if s.count(old)!=1: print(f"ERR [{tag}] count={s.count(old)}"); sys.exit(1)
    s=s.replace(old,new); print(f"ok {tag}")

# --- E1: extend S2 computation (KO + Nipponbare RBH copy) ---
old1='''# pathway copy number
steps={}
with open(PTSV) as fh:
    r=csv.reader(fh,delimiter='\\t'); next(r)
    for row in r:
        if not row or not row[0].strip(): continue
        steps.setdefault(row[0].strip(),0); steps[row[0].strip()]+=1
path_rows=[[k.split()[0], v] for k,v in steps.items()]; total_genes=sum(v for _,v in path_rows)'''
new1='''# pathway copy number (+ KO + Nipponbare RBH orthologue count)
RBH=ROOT+"/07_分析_Os02g基因鉴定/18_通路拷贝数_A11/best.tsv"
zn2nip={}
try:
    for ln in open(RBH):
        f=ln.rstrip("\\n").split("\\t")
        if len(f)>=2: zn2nip[f[0].split(".")[0]]=f[1]
except FileNotFoundError: pass
step_genes={}; step_ko={}
with open(PTSV) as fh:
    r=csv.reader(fh,delimiter='\\t'); next(r)
    for row in r:
        if not row or not row[0].strip(): continue
        st=row[0].strip(); step_genes.setdefault(st,[]).append(row[2].strip() if len(row)>2 else "")
        if len(row)>1 and row[1].strip(): step_ko.setdefault(st,row[1].strip())
path_rows=[]
for st,genes in step_genes.items():
    nip=len({zn2nip[g] for g in genes if g in zn2nip})
    path_rows.append([st.split()[0], step_ko.get(st,""), len(genes), nip if nip else "n.d."])
total_genes=sum(r[2] for r in path_rows)
total_nip=len({zn2nip[g] for gs in step_genes.values() for g in gs if g in zn2nip})'''
rep(old1,new1,"E1_S2compute")

# --- E2a: S2 docx table (new columns) ---
rep('''H(d,"Supplementary Table S2. Flavonoid/anthocyanin pathway gene copy number (KEGG orthology)",11)
tbl(d,["Pathway step","Copy number"],path_rows+[["TOTAL",total_genes]],fs=8.5)''',
'''H(d,"Supplementary Table S2. Flavonoid/anthocyanin pathway gene copy number: ZN65 vs Nipponbare (KEGG orthology)",11)
tbl(d,["Pathway step","KEGG KO","ZN65 copy","Nipponbare copy (RBH)"],path_rows+[["TOTAL","",total_genes,total_nip]],fs=8.5)
P(d,"Nipponbare copy number is the count of distinct IRGSP-1.0 reciprocal-best-hit orthologues (BLASTp) per step; "
  "a full OrthoFinder orthogroup analysis is the planned refinement. Copy numbers are from KEGG-orthology-annotated "
  "models and may include paralogues or pseudogenes; high-confidence single-copy curation is pending. ZN65 and "
  "Nipponbare have closely matched complements, indicating the pathway is structurally conserved rather than "
  "expanded.",size=9,color=GREY)''',"E2a_S2docx")

# --- E2b: S2 xlsx sheet ---
rep('''ws=wb.create_sheet("S2_pathway_copy_number"); ws.append(["Pathway step","Copy number"])
for r in path_rows: ws.append(r)
ws.append(["TOTAL",total_genes])''',
'''ws=wb.create_sheet("S2_pathway_copy_number"); ws.append(["Pathway step","KEGG KO","ZN65 copy","Nipponbare copy (RBH)"])
for r in path_rows: ws.append(r)
ws.append(["TOTAL","",total_genes,total_nip])''',"E2b_S2xlsx")

# --- E3: S1 footnote (docx) ---
rep('''tbl(d,["Variation type","Count","Length ref (bp)","Length qry (bp)"],sv_rows,fs=8.5)''',
'''tbl(d,["Variation type","Count","Length ref (bp)","Length qry (bp)"],sv_rows,fs=8.5)
P(d,"Counts are per-event-type SyRI annotations (minimap2 asm5). SyRI tabulates translocations (TRANS) and "
  "inverted translocations (INVTR), and duplications (DUP) and inverted duplications (INVDP), separately; the "
  "value quoted in the main text as 449 translocations is the aggregate TRANS+INVTR (198+251). SNP and small "
  "indels are sequence-level; inversions/translocations/duplications are structural. Lengths are in base pairs on "
  "the reference (Nipponbare) and query (ZN65) coordinates. Robustness to alignment stringency (asm5/10/20) and an "
  "orthogonal Assemblytics call are reported in the main-text Methods.",size=9,color=GREY)''',"E3_S1note")

# --- E4: S3 header + note strengthen ---
rep('''H(d,"Supplementary Table S3. Pigmentation-locus & MBW gene expression (vegetative tissue; raw counts and TPM)",11)''',
'''H(d,"Supplementary Table S3. Pigmentation-locus & MBW gene expression in VEGETATIVE (non-pericarp) tissue \\u2014 raw counts and TPM",11)''',"E4a_S3hdr")
rep('''P(d,"Note: vegetative-stage tissue (SD seedling, TG tillering); OsDFR and Rc are silent and OsTTG1 constitutive, "
  "consistent with pericarp-restricted pigmentation. Pericarp-stage data are the gating dataset G2.",size=9,color=GREY)''',
'''P(d,"IMPORTANT: tissue is VEGETATIVE and NON-PERICARP \\u2014 seedling (SD, n=3) and tillering (TG, n=2; TG2 library "
  "absent). No biological replicate exists for some contrasts, so no inferential statistics are reported here; values "
  "are raw counts and TPM only. OsDFR and Rc are silent and OsTTG1 is constitutive, consistent with "
  "pericarp-restricted pigmentation, but this baseline cannot demonstrate pericarp regulation. Pericarp "
  "developmental-stage expression is the gating dataset (G2). HISAT2 mapping 96.3\\u201397.8%.",size=9,color=GREY)''',"E4b_S3note")

# --- E5: S4 note strengthen ---
rep('''P(d,"Rc loss-of-function lesion in white rice (14-bp exon-6 deletion) per Sweeney et al. 2006. Nipponbare OsB2/OsDFR "
  "models may be truncated; Rc is the cleanest evidence line.",size=9,color=GREY)''',
'''P(d,"\\u201coverlap identity\\u201d is the percent identity over the aligned (overlapping) protein region only (e.g. Rc "
  "98.9% over the shared N-terminus), NOT full-length identity. CAVEAT: the Nipponbare OsB2 (180 aa) and OsDFR (284 aa) "
  "gene models may themselves be truncated/mis-annotated, so the \\u201cZN65 full-length vs NIP short\\u201d calls for "
  "OsB2/OsDFR are model-dependent; the cleanest, model-independent evidence is Rc, where white-pericarp Nipponbare "
  "carries the documented 14-bp exon-6 deletion removing the C-terminal bHLH (Sweeney et al. 2006) and ZN65 a "
  "full-length functional Rc.",size=9,color=GREY)''',"E5_S4note")

open(BP,"w",encoding="utf-8").write(s)
print("S1-S4 revisions applied to build_pkg.py")
