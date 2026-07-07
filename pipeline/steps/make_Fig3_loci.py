import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import pubstyle; font,C=pubstyle.apply()
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

GFF="/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.gff"
RM="/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/02.annotation/01.repeat/repeatmasker.gff"
# locus: (gene_id, chrom, win_start, win_end, title, note)
LOCI=[("ZN654G2687","Chr4",27995000,28028800,"Kala4 / OsB2  (Chr4, −)","TE-laden 14.4-kb intron (Gypsy) + ZN65-specific promoter insertion"),
      ("ZN657G0823","Chr7",6469000,6481500,"Rc  (Chr7, +)","Stowaway/MITE in promoter & introns; ~19% of locus non-aligned"),
      ("ZN651G2772","Chr1",27072000,27080500,"Kala1 / OsDFR  (Chr1, +)","Helitron cluster + Stowaway flanking the promoter"),
      ("ZN656G0716","Chr6",5189500,5195000,"OsC1 / Kala3  (Chr6, +)","colinear, 99.8% identity — TEs only in flanks (conserved)")]
def exons(gid):
    span=None; ex=[]
    for ln in open(GFF):
        if gid not in ln: continue
        f=ln.split("\t")
        if len(f)<8: continue
        if f[2]=="mRNA": span=(int(f[3]),int(f[4]),f[6])
        elif f[2] in ("CDS","exon"): ex.append((int(f[3]),int(f[4])))
    return span, sorted(set(ex))
def tes(chrom,a,b):
    out=[]
    for ln in open(RM):
        if not ln.startswith(chrom+"\t"): continue
        f=ln.split("\t")
        if len(f)<9: continue
        s,e=int(f[3]),int(f[4])
        if e<a or s>b: continue
        cls=""
        for kv in f[8].split(";"):
            if kv.startswith("Class="): cls=kv[6:]
        out.append((s,e,cls))
    return out
# TE class -> colour
def teclass(cls):
    c=cls.upper()
    if "GYPSY" in c: return ("LTR/Gypsy",C['vermilion'])
    if "COPIA" in c or c.startswith("LTR"): return ("LTR/Copia/other",C['orange'])
    if "LINE" in c: return ("LINE",'#8C4A2F')
    if "MUDR" in c or "MULE" in c: return ("DNA/MULE-MuDR",C['purple'])
    if "HELITRON" in c or c.startswith("RC"): return ("RC/Helitron",C['green'])
    if "STOWAWAY" in c or "TCMAR" in c or "MITE" in c: return ("DNA/Stowaway-MITE",C['sky'])
    if "SINE" in c: return ("SINE",C['rose'])
    return ("DNA/other",C['yellow'])

fig=plt.figure(figsize=(7.6,7.2))
gs=GridSpec(len(LOCI),1,hspace=0.95,figure=fig)
seen={}
for i,(gid,chrom,a,b,title,note) in enumerate(LOCI):
    ax=fig.add_subplot(gs[i]); ax.set_xlim(a/1000,b/1000); ax.set_ylim(-1.45,1.2); ax.axis('off')
    span,ex=exons(gid); gs0,gs1,strand=span
    sv=("OsC1" not in title)
    titlecol=C['vermilion'] if sv else C['blue']
    ax.text((a)/1000,0.95,title,fontsize=8.2,fontweight='bold',color=titlecol,ha='left')
    ax.text((a)/1000,0.62,note,fontsize=6.1,color='#555',ha='left')
    # gene model
    ax.add_line(Line2D([gs0/1000,gs1/1000],[0,0],color='#777',lw=1.0,zorder=2))
    for (es,ee) in ex:
        ax.add_patch(Rectangle((es/1000,-0.16),max((ee-es)/1000,0.03),0.32,facecolor=C['blue'],edgecolor='white',lw=0.3,zorder=4))
    # strand arrow
    if strand=='+':
        ax.annotate('',xy=(gs1/1000+0.15,0),xytext=(gs1/1000,0),arrowprops=dict(arrowstyle='-|>',color='#444',lw=1.0))
        ax.text(gs0/1000-0.05,0.30,"5′",fontsize=6,color='#333',ha='right')
    else:
        ax.annotate('',xy=(gs0/1000-0.15,0),xytext=(gs0/1000,0),arrowprops=dict(arrowstyle='-|>',color='#444',lw=1.0))
        ax.text(gs1/1000+0.05,0.30,"5′",fontsize=6,color='#333',ha='left')
    # TE track (below)
    yT=-0.82
    for s,e,cls in tes(chrom,a,b):
        lab,col=teclass(cls); seen[lab]=col
        ax.add_patch(Rectangle((s/1000,yT-0.16),max((e-s)/1000,0.04),0.32,facecolor=col,edgecolor='none',alpha=0.92,zorder=3))
    ax.text((a)/1000,yT+0.30,"transposable elements",fontsize=5.6,color='#888',ha='left',style='italic')
    # 1-kb scale bar
    ax.plot([b/1000-1.05,b/1000-0.05],[-1.20,-1.20],color='#333',lw=1.3)
    ax.text(b/1000-0.55,-1.30,"1 kb",fontsize=5.8,ha='center',va='top',color='#333')
    pubstyle.panel(ax,'abcd'[i],x=-0.02,y=1.02)
# legend (TE classes seen)
handles=[Rectangle((0,0),1,1,fc=C['blue'],label='gene exon (CDS)')]+[Rectangle((0,0),1,1,fc=c,label=l) for l,c in seen.items()]
fig.legend(handles=handles,loc='lower center',ncol=4,frameon=False,fontsize=6.2,bbox_to_anchor=(0.5,-0.02))
fig.suptitle("Gene-model and transposon architecture of the four pericarp-pigmentation loci (ZN65; per-panel scale differs, 1-kb bar shown)",
             fontsize=8.6,y=0.995,fontweight='bold')
plt.savefig("Fig3_loci.png",bbox_inches='tight'); plt.savefig("Fig3_loci.pdf",bbox_inches='tight')
print("saved Fig3_loci composite; TE classes:",list(seen))
