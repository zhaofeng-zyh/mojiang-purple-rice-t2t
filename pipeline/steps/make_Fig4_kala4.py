import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import pubstyle; font,C=pubstyle.apply()
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

GFF="/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.gff"
A10="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/A10/PA_table_ortho.tsv"
# --- auto-read OsB2 (ZN654G2687) exon model from GFF ---
mrna=None; exons=[]
for ln in open(GFF):
    if "ZN654G2687" not in ln: continue
    f=ln.split("\t")
    if len(f)<8: continue
    if f[2]=="mRNA": mrna=(int(f[3]),int(f[4]),f[6])
    elif f[2] in ("CDS","exon"): exons.append((int(f[3]),int(f[4])))
exons=sorted(set(exons)); g0,g1,strand=mrna
# verified RepeatMasker TE coords (结果_Kala4致色变异与定年.md)
gyp=(27999746,28012241); lin=(28014081,28020405)
prom=[(28021234,28021389,C['purple']),(28024034,28024144,C['green']),(28024370,28024447,C['sky']),
      (28025332,28025573,C['green']),(28025627,28025736,C['orange'])]
ins0,ins1=28022131,28028000   # A10 ZN65-specific insert (~5.9 kb; matches 'not aligned to Nipponbare')
x0,x1=27995500,28028800
mb=lambda v:v/1e6

fig=plt.figure(figsize=(8.6,7.2))
gs=GridSpec(2,1,height_ratios=[1.7,0.85],hspace=1.05,figure=fig)
ax=fig.add_subplot(gs[0]); ax.set_xlim(mb(x0),mb(x1)); ax.set_ylim(-3.1,2.55); ax.axis('off')
ax.axvspan(mb(g1),mb(x1),color='#FDF3E3',zorder=0)
ax.text(mb((g1+x1)/2),2.30,"5′ promoter / upstream  (− strand gene)",ha='center',fontsize=7.0,color='#9A6A00')
ax.text(mb((g0+g1)/2),2.30,"OsB2 / Kala4  (ZN654G2687, − strand, ~24 kb)",ha='center',fontsize=7.0,color=C['blue'])

# intron baseline + exon boxes (real gene model)
ax.add_line(Line2D([mb(g0),mb(g1)],[0,0],color='#777',lw=1.0,zorder=1))
for (es,ee) in exons:
    ax.add_patch(Rectangle((mb(es),-0.16),max(mb(ee-es),0.00025),0.32,facecolor=C['blue'],edgecolor='white',lw=0.3,zorder=4))
# big first intron (contains Gypsy): between the two 5' exon clusters
ex_sorted=sorted(exons); 
# find largest intron gap
gaps=[(ex_sorted[i+1][0]-ex_sorted[i][1], ex_sorted[i][1], ex_sorted[i+1][0]) for i in range(len(ex_sorted)-1)]
big=max(gaps); 
ax.annotate(f"largest intron (~{big[0]/1000:.1f} kb) — harbours the intronic Gypsy",
            xy=(mb((big[1]+big[2])/2),0.02),xytext=(mb((big[1]+big[2])/2),0.95),
            fontsize=6.0,color='#444',ha='center',arrowprops=dict(arrowstyle='-',color='#aaa',lw=0.5))
# ATG / direction (− strand: start at high-coord end)
ax.add_patch(FancyArrow(mb(g1),0,-0.0011,0,width=0,head_width=0.30,head_length=0.0006,length_includes_head=True,color='#444',zorder=5))
ax.text(mb(g1)+0.0002,0.42,"ATG / TSS (5′)",ha='left',va='bottom',fontsize=6.2,color='#333')
ax.text(mb(g0)-0.0002,0.42,"3′",ha='right',va='bottom',fontsize=6.2,color='#333')

# TE track (below gene)
yT=-1.05
ax.add_patch(Rectangle((mb(gyp[0]),yT-0.16),mb(gyp[1]-gyp[0]),0.32,facecolor=C['vermilion'],edgecolor='white',lw=0.4,zorder=3))
ax.add_patch(Rectangle((mb(lin[0]),yT-0.16),mb(lin[1]-lin[0]),0.32,facecolor=C['orange'],edgecolor='white',lw=0.4,zorder=3))
ax.text(mb((gyp[0]+gyp[1])/2),yT+0.46,"Gypsy LTR-RT (RETRO2B, ~12.5 kb, in intron)",ha='center',fontsize=6.2,color=C['vermilion'],fontweight='bold')
ax.text(mb((lin[0]+lin[1])/2),yT+0.46,"LINE-1 (~6.3 kb)",ha='center',fontsize=5.8,color=C['orange'])
for s,e,col in prom:
    ax.add_patch(Rectangle((mb(s),yT-0.16),max(mb(e-s),0.00022),0.32,facecolor=col,edgecolor='white',lw=0.3,zorder=3))
ax.text(mb((prom[0][0]+prom[-1][1])/2),yT-0.40,"promoter TE cluster (SINE·Helitron·hAT·PIF·LINE-1)",ha='center',va='top',fontsize=5.7,color='#555')
# dashed link: Gypsy sits in the big intron
ax.plot([mb(gyp[0]),mb(gyp[0])],[0,yT+0.16],color='#bbb',lw=0.5,ls=':'); ax.plot([mb(gyp[1]),mb(gyp[1])],[0,yT+0.16],color='#bbb',lw=0.5,ls=':')

# insert bracket (A10) + dating
ax.annotate('',xy=(mb(ins0),1.42),xytext=(mb(ins1),1.42),arrowprops=dict(arrowstyle='|-|',color=C['red'],lw=1.1))
ax.text(mb((ins0+ins1)/2),1.74,f"ZN65-specific promoter insertion (~{(ins1-ins0)/1000:.1f} kb)",ha='center',va='bottom',fontsize=5.8,color=C['red'],fontweight='bold')
ax.text(mb((gyp[0]+gyp[1])/2),-2.05,"intronic Gypsy LTR age ≈ 0.2 Myr (95% CI 0.04–0.39; K2P, mu=1.3e-8)",ha='center',va='center',fontsize=6.0,color=C['vermilion'],bbox=dict(boxstyle='round,pad=0.3',fc='white',ec=C['vermilion'],lw=0.6))
# scale axis
ax.plot([mb(x0),mb(x1)],[-2.45,-2.45],color='#333',lw=0.7)
for xt in [27.996,28.004,28.012,28.020,28.028]:
    ax.plot([xt,xt],[-2.45,-2.49],color='#333',lw=0.6); ax.text(xt,-2.52,f"{xt:.3f}",ha='center',va='top',fontsize=6,color='#333')
ax.text(mb((x0+x1)/2),-2.82,"ZN65 chromosome 4 position (Mb)",ha='center',va='top',fontsize=7.0)
pubstyle.panel(ax,'a',x=-0.02,y=1.10)

# ---- panel b: quantitative ortholog-anchored presence (A10) ----
axb=fig.add_subplot(gs[1])
rows=[]
for ln in open(A10):
    f=ln.rstrip("\n").split("\t")
    if f[0]=="genome" or len(f)<6: continue
    rows.append((f[0],f[2],float(f[4])))  # genome, pericarp, insert_localcov%
order=["ZN65","Nipponbare","MH63","ZS97","CempoIreng","N22","Orufipogon"]
rows=sorted(rows,key=lambda r:order.index(r[0]) if r[0] in order else 99)
labs=[f"{g} ({p})" for g,p,_ in rows]; vals=[v for _,_,v in rows]
import numpy as np
y=np.arange(len(rows))[::-1]
cols=[C['red'] if g=="ZN65" else C['grey'] for g,_,_ in rows]
axb.barh(y,vals,color=cols,height=0.6,edgecolor='white',linewidth=0.3)
for yi,v in zip(y,vals): axb.text(v+2,yi,("present" if v>=60 else "absent")+f" ({v:.0f}%)",va='center',fontsize=6.0,color=('#B0271A' if v>=60 else '#777'))
axb.set_yticks(y); axb.set_yticklabels(labs,fontsize=6.4)
axb.set_xlim(0,118); axb.set_xlabel("ZN65-specific insertion: local coverage at the ORTHOLOGOUS OsB2 locus (%)",fontsize=6.6)
axb.set_title("Insertion is ZN65-lineage-specific — absent from the orthologous locus in indica, aus, wild and black rice",fontsize=6.8,loc='left',pad=4)
axb.text(0.5,-1.15,"Method: 34-kb conserved OsB2 region anchors the orthologous locus per genome; the 5.9-kb insert is scored only within that interval "
         "(avoids dispersed-TE false positives). Caveat: one genome per accession/species — intraspecific polymorphism not assessed; "
         "population/pan-genome sampling is the ideal extension.",transform=axb.transAxes,fontsize=5.4,color='#666',ha='left',va='top',wrap=True)
for s in ('top','right'): axb.spines[s].set_visible(False)
pubstyle.panel(axb,'b',x=-0.02,y=1.12)

fig.suptitle("Retrotransposon architecture and lineage-specific origin of the Kala4/OsB2 pigmentation allele in ZN65",fontsize=9.0,y=1.0,fontweight='bold')
plt.savefig("Fig4_kala4_architecture.png",bbox_inches='tight'); plt.savefig("Fig4_kala4_architecture.pdf",bbox_inches='tight')
print(f"saved Fig4: OsB2 {len(exons)} exons, big intron {big[0]/1000:.1f}kb, insert {(ins1-ins0)/1000:.1f}kb")
