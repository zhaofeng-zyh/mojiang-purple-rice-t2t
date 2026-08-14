# ⛔ SUPERSEDED —— 本脚本按**已被推翻的口径**作图，入口硬失败（2026-08-01 加装）。
# 第 29 行 `specific=73.3` **参与算术**：nonte = 73.3 − 45.12 = 28.18；
#   按现行口径 **73.2 Mb**（= NOTAL 70,339,400 + INS 2,890,667 = 73,230,067 bp）
#   应为 28.08 —— 非 TE 段差 0.10 Mb。
# **投稿用的 Figure 2C 不是本脚本产出**：现行生成器是
#   `12_论文Paper1_Manuscript/_投稿包_TheCropJournal_20260724/_figures_src/make_Fig2_SV_v2.py`，
#   它用 73.2（`total = sum(seg_val)`），投稿图与稿件口径一致。本文件保留仅作留痕。
import os as _os, sys as _sys
if _os.environ.get("ALLOW_SUPERSEDED") != "1":
    _sys.exit("⛔ 本脚本用已被推翻的 73.3 Mb 作图（现行 73.2 Mb），"
              "现行生成器见 _投稿包_TheCropJournal_20260724/_figures_src/make_Fig2_SV_v2.py。"
              "确需留痕重现请设 ALLOW_SUPERSEDED=1。")
import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import pubstyle; font,C=pubstyle.apply()
import numpy as np, matplotlib.pyplot as plt, matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

R="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP/"
HERE=os.path.dirname(os.path.abspath(__file__))
sv={}
for ln in open(R+"zn65nip_syri.summary"):
    if ln.startswith("#") or not ln.strip(): continue
    p=ln.rstrip().split("\t")
    try: sv[p[0]]=int(p[1])
    except: pass
# structural-variant length distribution (block sizes on ref, from syri_structural.out)
svlen=[]
struct={"INV","TRANS","DUP","INVDP","INVTR"}
try:
    for ln in open(R+"syri_structural.out"):
        f=ln.rstrip("\n").split("\t")
        if len(f)<11 or f[10] not in struct: continue
        try:
            a,b=int(f[1]),int(f[2]); L=abs(b-a)+1
            if L>=1: svlen.append(L)
        except: pass
except FileNotFoundError: pass
# TE composition of ZN65-specific 73.3 Mb (A16); add non-TE contrast
te=[("LTR/Gypsy",26.91),("DNA transposons",11.11),("LTR/Copia",3.00),("LINE",1.33),("LTR/other",0.63),("SINE",0.10)]
te_total=45.12; specific=73.3; nonte=specific-te_total

fig=plt.figure(figsize=(7.8,8.6))
gs=GridSpec(2,3,height_ratios=[2.15,1.0],hspace=0.26,wspace=0.52,figure=fig)

# (a) genome-wide synteny ribbon (plotsr) as MAIN panel
axa=fig.add_subplot(gs[0,:])
img=mpimg.imread(os.path.join(HERE,"synteny_ribbon.png"))
axa.imshow(img,aspect='auto'); axa.axis('off')
axa.set_title("ZN65 vs Nipponbare whole-genome synteny (SyRI/plotsr): grey, syntenic; "
              "orange, inversion (note Chr6); green, translocation; blue, duplication",
              fontsize=7.2,loc='left',pad=4)
pubstyle.panel(axa,'a',x=-0.02,y=1.01)

# (b) variant spectrum, split: sequence-level vs structural
axb=fig.add_subplot(gs[1,0])
seq_order=[('SNPs',C['blue']),('Insertions',C['sky']),('Deletions',C['sky'])]
str_order=[('Highly diverged',C['vermilion']),('Inversions',C['orange']),('Translocations',C['vermilion']),
           ('Duplications (query)',C['green']),('Duplications (reference)',C['green'])]
short={'Duplications (query)':'Dup.(qry)','Duplications (reference)':'Dup.(ref)','Highly diverged':'Highly div.','Translocations':'Transloc.','Inversions':'Inversion'}
rows=[(k,c,'seq') for k,c in seq_order if k in sv]+[(k,c,'str') for k,c in str_order if k in sv]
y=np.arange(len(rows))[::-1]
axb.barh(y,[sv[k] for k,_,_ in rows],color=[c for _,c,_ in rows],height=0.66,edgecolor='white',linewidth=0.3)
for yi,(k,_,_) in zip(y,rows): axb.text(sv[k]*1.3,yi,f"{sv[k]:,}",va='center',fontsize=5.4)
axb.set_xscale('log'); axb.set_xlim(1,max(sv.values())*6)
axb.set_yticks(y); axb.set_yticklabels([short.get(k,k) for k,_,_ in rows],fontsize=6.2)
# divider between sequence-level and structural groups
nseq=sum(1 for r in rows if r[2]=='seq')
axb.axhline(len(rows)-nseq-0.5,color='#bbb',lw=0.7,ls='--')
axb.text(1.4,len(rows)-0.5,"sequence-level",fontsize=5.6,color=C['blue'],va='center',style='italic')
axb.text(1.4,len(rows)-nseq-1.0,"structural",fontsize=5.6,color=C['vermilion'],va='top',style='italic')
axb.set_xlabel("Count (log)"); axb.set_title("Variant spectrum",fontsize=7.6,loc='left',pad=4)
pubstyle.panel(axb,'b',x=-0.34)
for s in ('top','right'): axb.spines[s].set_visible(False)

# (c) TE composition + non-TE contrast of ZN65-specific sequence
axc=fig.add_subplot(gs[1,1])
labels=[l for l,_ in te]+["non-TE (unique)"]; mb=[v for _,v in te]+[nonte]
palette=[C['vermilion'],C['orange'],C['rose'],C['sky'],C['yellow'],C['green'],C['grey']]
yy=np.arange(len(labels))[::-1]
axc.barh(yy,mb,color=palette,height=0.66,edgecolor='white',linewidth=0.3)
for yi,v in zip(yy,mb): axc.text(v+0.4,yi,f"{v:.1f}",va='center',fontsize=5.8)
axc.set_yticks(yy); axc.set_yticklabels(labels,fontsize=6.2)
axc.set_xlabel("Sequence (Mb)"); axc.set_xlim(0,max(mb)*1.18)
axc.set_title(f"ZN65-specific seq ({specific:.0f} Mb)\nTE {te_total:.1f} Mb (62%) + non-TE {nonte:.1f} Mb",fontsize=6.9,loc='left',pad=4)
pubstyle.panel(axc,'c',x=-0.30)
for s in ('top','right'): axc.spines[s].set_visible(False)

# (d) structural-variant length distribution
axd=fig.add_subplot(gs[1,2])
if svlen:
    sl=np.array(svlen); sl=sl[sl>=10]
    bins=np.logspace(np.log10(max(sl.min(),10)),np.log10(sl.max()),22)
    axd.hist(sl,bins=bins,color=C['purple'],edgecolor='white',linewidth=0.3)
    axd.set_xscale('log'); axd.set_xlabel("SV length (bp, log)"); axd.set_ylabel("Count")
    axd.set_title(f"Structural-variant sizes\n(INV/TRANS/DUP; n={len(sl):,})",fontsize=6.9,loc='left',pad=4)
else:
    axd.text(0.5,0.5,"SV length data\nunavailable",ha='center',va='center',fontsize=7); axd.axis('off')
pubstyle.panel(axd,'d',x=-0.30)
for s in ('top','right'): axd.spines[s].set_visible(False)

fig.suptitle("Genome-wide structural divergence and transposon-driven expansion of ZN65 vs Nipponbare",
             fontsize=9.6,y=0.995,fontweight='bold')
plt.savefig("Fig2_genomewide_SV.png",bbox_inches='tight'); plt.savefig("Fig2_genomewide_SV.pdf",bbox_inches='tight')
print("saved Fig2 (synteny main + split spectrum + TE/nonTE + SV-length); n_svlen=%d"%len(svlen))
