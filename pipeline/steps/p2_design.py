# ⛔ SUPERSEDED —— 本脚本会按**已被推翻的口径**再生产物，入口硬失败（2026-07-29 加装）
#
# 已被推翻的旧口径：「OsB2 启动子有一段 ~5.9 kb / ~6.3 kb 的 ZN65 特异插入，白等位缺失，
# 可用存在-缺失共显性分型」。现行口径：OsB2 **上游**是 ~5.87 kb block A 的 LINE 关联
# **复合局部片段重复** —— ZN65 3 拷贝（A1/A2/A3），日本晴/MH63/ZS97/N22/Cempo Ireng/
# 野生 O. rufipogon **6 个对照各单拷贝**（99.3–100%）；该片段**在所有 7 个受检基因组中都存在**，
# 差异只是拷贝数 3 vs 1，**不存在「对照缺失」**，故**不能做存在-缺失共显性分型**。
# 接头缺 TSD 与完整末端重复 → **不是逆转座子插入**。内含子 LINE1-11_OS 实测 **5,934 bp**（非 6.3 kb）；
# 三拷贝区 A1 起–A3 止 **115,450 bp**（非 ~119 kb）。
#
# 保留本文件仅为留痕。若确需运行（例如复现旧图以作对照），显式设
# `ALLOW_SUPERSEDED=1` —— 但产物**不得**作为现行结论使用。
import os as _os, sys as _sys
if _os.environ.get("ALLOW_SUPERSEDED") != "1":
    _sys.exit("⛔ 本脚本按已被推翻的口径产出（详见文件头注释）。"
              "现行口径见 12_论文Paper1_Manuscript/_投稿包_TheCropJournal_20260724/。"
              "确需留痕重现请设 ALLOW_SUPERSEDED=1。")
import subprocess, primer3, re
B="/mnt/data2/墨江紫米研究"
REF=B+"/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa"
W=B+"/13_Paper2_机制研究方案"
def faidx(region):
    out=subprocess.run(["seqkit","faidx",REF,region],capture_output=True,text=True).stdout
    return "".join(l.strip() for l in out.split("\n")[1:]).upper()
def rc(s):
    return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def gc(s): return 100*(s.count("G")+s.count("C"))/len(s)
# blast db (genome) for uniqueness
subprocess.run(["makeblastdb","-in",REF,"-dbtype","nucl","-out","/tmp/gdb"],capture_output=True)
def n_hits(seq,ident=90,minlen=None):
    open("/tmp/q.fa","w").write(f">q\n{seq}\n")
    o=subprocess.run(["blastn","-task","blastn-short","-query","/tmp/q.fa","-db","/tmp/gdb",
        "-outfmt","6 length pident","-evalue","10","-perc_identity",str(ident)],capture_output=True,text=True).stdout
    ml=minlen or int(0.9*len(seq))
    return sum(1 for l in o.strip().split("\n") if l and int(l.split("\t")[0])>=ml)

OUT=[]
OUT.append("# Paper2 分子试剂设计（ZN65 T2T 基因组）—— 自动设计 + 特异性已查\n")

# ---------- 1. CRISPR sgRNA: OsB2 ZN654G2687, target exon1+exon2 (5'编码区) ----------
OUT.append("## 1. CRISPR-Cas9 sgRNA（OsB2 / ZN654G2687，靶 5' 编码区致移码）")
exons={"exon1":(28020506,28020642),"exon2":(28013739,28013993)}
def find_sgRNAs(chrom,s,e,label):
    region=faidx(f"{chrom}:{s-3}-{e+3}")  # pad for PAM
    cands=[]
    # sense: protospacer(20) + NGG
    for m in re.finditer(r'(?=([ACGT]{21}GG))',region):
        full=m.group(1); proto=full[:20]
        cands.append((proto,"+",full[20:23]))
    for m in re.finditer(r'(?=(CC[ACGT]{21}))',region):
        full=m.group(1); proto=rc(full[3:])
        cands.append((proto,"-","CCN->"+rc(full[:3])))
    res=[]
    for proto,strand,pam in cands:
        if "TTTT" in proto: continue
        g=gc(proto)
        if not (40<=g<=75): continue
        h=n_hits(proto,ident=95)
        res.append((proto,strand,round(g),h))
    # prefer unique (h==1), pick top2 by GC near 55
    res=sorted(set(res),key=lambda r:(r[3]!=1, abs(r[2]-55)))
    return res[:3]
for lab,(s,e) in exons.items():
    rs=find_sgRNAs("Chr4",s,e,lab)
    for i,(proto,strand,g,h) in enumerate(rs[:2],1):
        flag="唯一" if h==1 else f"⚠{h}处匹配(查脱靶)"
        OUT.append(f"- sgRNA-OsB2-{lab}-{i}: 5'-{proto}-3'  (链{strand}, GC{g}%, {flag}; 接头加 BsaI/BsmBI 按 pYLCRISPR/Cas9 Ma2015)")
OUT.append("> 设计参数: 20 nt protospacer + NGG PAM; GC 40–75%; 去 polyT(TTTT); 已 blast ZN65 基因组查唯一性。**订购前再用 CRISPR-P/CRISPOR 复核脱靶, 并确认靶区在 ZN65 无 SNP。**\n")

# ---------- 2. 插入断点/连接 基因型分型引物 (共分离 + 种质panel + KASP底物) ----------
OUT.append("## 2. OsB2 启动子插入 断点分型引物（共分离/种质panel，存在/缺失共显性）")
ins0,ins1=28022131,28028000
upF=faidx(f"Chr4:{ins0-350}-{ins0-1}")      # 上游单拷贝(插入5'外侧)
dnR=faidx(f"Chr4:{ins1+1}-{ins1+350}")       # 下游单拷贝(插入3'外侧)
ins5=faidx(f"Chr4:{ins0}-{ins0+350}")        # 插入体5'端内部
ins3=faidx(f"Chr4:{ins1-350}-{ins1}")        # 插入体3'端内部
def p3(seq,tag,want_left=True,want_right=True):
    g={'PRIMER_OPT_SIZE':22,'PRIMER_MIN_SIZE':19,'PRIMER_MAX_SIZE':25,'PRIMER_OPT_TM':60,'PRIMER_MIN_TM':58,'PRIMER_MAX_TM':63,
       'PRIMER_MIN_GC':40,'PRIMER_MAX_GC':60,'PRIMER_PRODUCT_SIZE_RANGE':[[100,320]],'PRIMER_NUM_RETURN':3,
       'PRIMER_PICK_LEFT_PRIMER':int(want_left),'PRIMER_PICK_RIGHT_PRIMER':int(want_right),'PRIMER_PICK_INTERNAL_OLIGO':0}
    r=primer3.bindings.design_primers({'SEQUENCE_ID':tag,'SEQUENCE_TEMPLATE':seq},g)
    return r
# single best primers from each window
def best_oligo(seq,side):
    g={'PRIMER_OPT_SIZE':22,'PRIMER_MIN_SIZE':19,'PRIMER_MAX_SIZE':25,'PRIMER_OPT_TM':60,'PRIMER_MIN_TM':58,'PRIMER_MAX_TM':63,
       'PRIMER_MIN_GC':40,'PRIMER_MAX_GC':62,'PRIMER_NUM_RETURN':1,'PRIMER_PICK_INTERNAL_OLIGO':0,
       'PRIMER_PICK_LEFT_PRIMER':int(side=='L'),'PRIMER_PICK_RIGHT_PRIMER':int(side=='R')}
    r=primer3.bindings.design_primers({'SEQUENCE_ID':'x','SEQUENCE_TEMPLATE':seq},g)
    k='PRIMER_LEFT_0_SEQUENCE' if side=='L' else 'PRIMER_RIGHT_0_SEQUENCE'
    return r.get(k,"NA")
Fout=best_oligo(upF,'L')      # 上游正向 (插入外侧5')
Rflank=best_oligo(dnR,'R')    # 下游反向 (插入外侧3')
Rins=best_oligo(ins5,'R')     # 插入内部反向 (5'连接诊断)
Fins=best_oligo(ins3,'L')     # 插入内部正向 (3'连接诊断)
OUT.append(f"- F-flank5 (插入上游, 正向): 5'-{Fout}-3'")
OUT.append(f"- R-flank3 (插入下游, 反向): 5'-{Rflank}-3'")
OUT.append(f"- R-ins5  (插入体5'端, 反向): 5'-{Rins}-3'")
OUT.append(f"- F-ins3  (插入体3'端, 正向): 5'-{Fins}-3'")
OUT.append("**判读**: ①F-flank5 + R-ins5 出带 = 含插入(紫等位 5'连接); ②F-ins3 + R-flank3 出带 = 含插入(3'连接); ③F-flank5 + R-flank3 仅在**无插入**(白等位, ~"+str(round((ins1-ins0)/1000,1))+"kb 太长不扩)出小带——共显性区分纯合紫/白与杂合。可据此做 KASP。\n")

# ---------- 3. MBW CDS 克隆引物 (Y2H/BiFC/原核表达) ----------
OUT.append("## 3. MBW 基因 CDS 克隆引物（Y2H/BiFC/Split-LUC/原核表达；建议 Gateway 或 同源重组）")
cds={}
cur=None
for ln in open(W+"/qPCR_target_CDS_for_primer_design.fasta"):
    if ln[0]==">": cur=ln[1:].split()[0].split(".")[0]; cds[cur]=[]
    else: cds[cur].append(ln.strip())
cds={k:"".join(v).upper() for k,v in cds.items()}
mbw={"OsB2(ZN654G2687)":"ZN654G2687","OsC1(ZN656G0716)":"ZN656G0716","OsTTG1(ZN652G3195)":"ZN652G3195"}
for name,gid in mbw.items():
    s=cds[gid]
    f=s[:22]; r=rc(s[-22:])           # ORF 两端 (含/不含终止子按融合方向)
    OUT.append(f"- {name}: F(ATG起)=5'-{f}-3'  R(含/去终止子)=5'-{r}-3'  (CDS {len(s)} bp; 加 attB1/attB2 或载体同源臂)")
OUT.append("> N端融合(BiFC-N/AD)用去终止子反向引物; 全长 ORF 已在 qPCR_target_CDS_*.fasta。\n")

# ---------- 4. EMSA / Y1H 顺式元件探针 (DFR=ZN651G2772, ANS=ZN652G3429 启动子) ----------
OUT.append("## 4. EMSA/Y1H 顺式元件探针（结构基因启动子中的 MYB/bHLH 结合位点）")
targets={"OsDFR(ZN651G2772)":("Chr1",15900232,15901825,"-"),"OsANS(ZN652G3429)":("Chr2",31708527,31711284,"+")}
# 取 TSS 上游 ~1.5kb 启动子
motifs={"G-box(bHLH)":"CACGTG","E-box(bHLH)":r"CA[ACGT][ACGT]TG","MRE/MBSI(MYB)":r"[ACT][ACG]C[CT]A[AC]C","AC-element(MYB)":r"ACC[AT]A[AC][CT]"}
for name,(chrom,gs,ge,strand) in targets.items():
    if strand=="+": pstart,pend=gs-1500,gs-1
    else: pstart,pend=ge+1,ge+1500
    prom=faidx(f"{chrom}:{pstart}-{pend}")
    if strand=="-": prom=rc(prom)   # 取转录方向
    OUT.append(f"### {name} 启动子(TSS上游~1.5kb)")
    found=[]
    for mlab,pat in motifs.items():
        for m in re.finditer(pat,prom):
            pos=m.start(); 
            if any(abs(pos-p)<6 for _,p in found): continue
            found.append((mlab,pos))
    # 取前若干个, 输出 ~46bp 探针(元件居中) + 突变对照
    for mlab,pos in found[:4]:
        a=max(0,pos-20); b=min(len(prom),pos+26)
        probe=prom[a:b]
        OUT.append(f"- {mlab} @+{pos}: 探针 5'-{probe}-3' (生物素/Cy5 标记双链; 突变对照=核心元件碱基替换; 冷探针竞争)")
OUT.append("> 探针为双链退火寡核苷酸, 5'生物素或 Cy5 标记; 每元件配**突变探针**(破坏核心序列)与**冷探针竞争**对照; bHLH(OsB2)需与 MYB(OsC1)共孵或仅测组合。\n")

open(W+"/Paper2_分子试剂设计清单.md","w").write("\n".join(OUT))
print("\n".join(OUT))
print("\n>>> saved Paper2_分子试剂设计清单.md")
