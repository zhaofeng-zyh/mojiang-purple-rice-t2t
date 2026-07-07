import primer3, subprocess, os, re
B="/mnt/data2/墨江紫米研究"
CDS_ALL=B+"/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.cds"
GFF=B+"/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.gff"
W=B+"/13_Paper2_机制研究方案"
genes=[("OsB2/Kala4","ZN654G2687","bHLH 调控(核心)"),("OsC1/Kala3","ZN656G0716","R2R3-MYB 调控"),
       ("OsTTG1","ZN652G3195","WD40 (MBW)"),("OsDFR/Kala1","ZN651G2772","结构 committed(单拷贝)"),
       ("ANS/LDOX","ZN652G3429","结构 committed(多拷贝代表)"),("UFGT/3GT","ZN656G0641","结构 committed(多拷贝代表)"),
       ("Actin(内参)","ZN653G3950","reference"),("UBQ(内参)","ZN653G1043","reference"),("eEF1a(内参)","ZN653G0605","reference")]
# CDS seqs
seqs={}; cur=None
for ln in open(W+"/qPCR_target_CDS_for_primer_design.fasta"):
    if ln[0]==">": cur=ln[1:].split()[0].split(".")[0]; seqs[cur]=[]
    else: seqs[cur].append(ln.strip())
seqs={k:"".join(v) for k,v in seqs.items()}
# exon junctions in CDS coords (strand-aware)
def junctions(gid):
    ex=[]; strand="+"
    for ln in open(GFF):
        if gid not in ln: continue
        f=ln.split("\t")
        if len(f)<8: continue
        if f[2]=="CDS": ex.append((int(f[3]),int(f[4]))); strand=f[6]
    ex=sorted(set(ex), reverse=(strand=="-"))
    js=[]; c=0
    for s,e in ex[:-1]: c+=(e-s+1); js.append(c)
    return js, len(ex)
GA={'PRIMER_PICK_LEFT_PRIMER':1,'PRIMER_PICK_RIGHT_PRIMER':1,'PRIMER_PICK_INTERNAL_OLIGO':0,
 'PRIMER_OPT_SIZE':20,'PRIMER_MIN_SIZE':18,'PRIMER_MAX_SIZE':22,
 'PRIMER_OPT_TM':60.0,'PRIMER_MIN_TM':58.0,'PRIMER_MAX_TM':62.0,'PRIMER_PAIR_MAX_DIFF_TM':2.0,
 'PRIMER_MIN_GC':40.0,'PRIMER_MAX_GC':60.0,'PRIMER_GC_CLAMP':1,'PRIMER_MAX_POLY_X':4,
 'PRIMER_PRODUCT_SIZE_RANGE':[[90,160]],'PRIMER_NUM_RETURN':8}
def specific(fwd,rev,gid):
    open("/tmp/q.fa","w").write(f">F\n{fwd}\n>R\n{rev}\n")
    out=subprocess.run(["blastn","-task","blastn-short","-query","/tmp/q.fa","-db","/tmp/cdsdb",
        "-outfmt","6 qseqid sseqid pident length","-evalue","1","-perc_identity","90"],
        capture_output=True,text=True).stdout
    genes_hit=set()
    for l in out.strip().split("\n"):
        if not l: continue
        p=l.split("\t")
        if int(p[3])>=int(0.9*len(fwd)): genes_hit.add(p[1].split(".")[0])
    return len(genes_hit)<=1, len(genes_hit)
# build cds blast db
subprocess.run(["makeblastdb","-in",CDS_ALL,"-dbtype","nucl","-out","/tmp/cdsdb"],capture_output=True)
rows=[]
for name,gid,role in genes:
    seq=seqs.get(gid,"")
    js,nex=junctions(gid)
    if not seq: rows.append((name,gid,"NO CDS","","","","","")); continue
    res=primer3.bindings.design_primers({'SEQUENCE_ID':gid,'SEQUENCE_TEMPLATE':seq},GA)
    n=res.get('PRIMER_PAIR_NUM_RETURNED',0)
    best=None
    for i in range(n):
        ls=res[f'PRIMER_LEFT_{i}'][0]; ll=res[f'PRIMER_LEFT_{i}'][1]
        rs=res[f'PRIMER_RIGHT_{i}'][0]  # right start (3' on template)
        amp_lo=ls; amp_hi=rs
        spans=any(amp_lo< j <amp_hi for j in js)   # amplicon crosses an exon junction (→ spans intron in gDNA)
        fwd=res[f'PRIMER_LEFT_{i}_SEQUENCE']; rev=res[f'PRIMER_RIGHT_{i}_SEQUENCE']
        ok,nhit=specific(fwd,rev,gid)
        cand=(spans,ok,i,fwd,rev,res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
              round(res[f'PRIMER_LEFT_{i}_TM'],1),round(res[f'PRIMER_RIGHT_{i}_TM'],1),
              round(res[f'PRIMER_LEFT_{i}_GC_PERCENT']),round(res[f'PRIMER_RIGHT_{i}_GC_PERCENT']),nhit)
        # prefer: specific + junction-spanning
        score=(ok*2)+(spans*1)
        if best is None or score>best[0]: best=(score,cand)
    if best is None: rows.append((name,gid,"no primer","","","","","")); continue
    _,(spans,ok,i,fwd,rev,psize,tmL,tmR,gcL,gcR,nhit)=best
    flag=("跨外显子" if spans else ("单外显子-需DNase" if nex==1 else "同外显子-需DNase")) + ("" if ok else f" ⚠多匹配{nhit}基因")
    rows.append((name,gid,fwd,rev,f"{psize}",f"{tmL}/{tmR}",f"{gcL}/{gcR}",flag))
# print table
print("基因\tZN65 ID\t上游引物F(5'->3')\t下游引物R(5'->3')\t产物bp\tTm F/R\tGC% F/R\t标注")
for r in rows: print("\t".join(map(str,r)))
# save TSV
open(W+"/qPCR_primers_synthesis_list.tsv","w").write(
  "Gene\tZN65_ID\tForward(5-3)\tReverse(5-3)\tProduct_bp\tTm_F/R\tGC_F/R\tNote\n"+
  "\n".join("\t".join(map(str,r)) for r in rows)+"\n")
print("\nsaved qPCR_primers_synthesis_list.tsv")
