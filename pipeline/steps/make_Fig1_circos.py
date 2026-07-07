import sys,os,math
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import pubstyle; font,C=pubstyle.apply()
import numpy as np, matplotlib.pyplot as plt, matplotlib as mpl
from matplotlib import colormaps
from matplotlib.colors import Normalize

D="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/11_基因组景观/"
chrlen={}
for ln in open(D+"genome.txt"):
    c,l=ln.split(); chrlen[c]=int(l)
chroms=[f"Chr{i}" for i in range(1,13)]
def load(fn,col):
    d={c:[] for c in chroms}
    for ln in open(D+fn):
        p=ln.split()
        if p[0] in d:
            try: d[p[0]].append((int(p[1]),float(p[col])))
            except: pass
    for c in d: d[c].sort()
    return d
gene=load("gene_density.txt",3); te=load("te_cov.txt",6); gc=load("gc.txt",3)
WIN=500000; total=sum(chrlen.values())
# pigmentation loci (genome-absolute), for outer-ring markers
LOCI=[("Chr4",27996908,"OsB2/Kala4"),("Chr1",27075322,"OsDFR/Kala1"),
      ("Chr6",5191615,"OsC1/Kala3"),("Chr7",6471567,"Rc")]
# GLOBAL per-track 2-98 percentile limits (comparable across chromosomes)
def glob_lim(track):
    allv=[v for c in chroms for _,v in track[c]]
    return np.percentile(allv,2),np.percentile(allv,98)
g_lo,g_hi=glob_lim(gene); t_lo,t_hi=glob_lim(te); gc_lo,gc_hi=glob_lim(gc)
def gnorm(v,lo,hi): return float(np.clip((v-lo)/((hi-lo) or 1e-9),0,1))
# approximate centromere per chr = window of maximum TE coverage
cen={}
for c in chroms:
    if te[c]: cen[c]=max(te[c],key=lambda x:x[1])[0]+WIN/2

GAP=math.radians(2.2); usable=2*math.pi-len(chroms)*GAP
spans={}; cur=math.pi/2
for c in chroms:
    ext=usable*chrlen[c]/total; spans[c]=(cur,cur-ext); cur=cur-ext-GAP
def th(c,pos):
    s,e=spans[c]; return s+(e-s)*(pos/chrlen[c])

fig=plt.figure(figsize=(8.2,9.0)); ax=fig.add_subplot(111,projection='polar')
ax.set_theta_zero_location('N'); ax.set_theta_direction(-1); ax.axis('off'); ax.set_ylim(0,1)
R_IDEO=(0.90,0.965); R_GENE=(0.755,0.875); R_TE=(0.61,0.73); R_GC=(0.45,0.58)
def ring_heat(track,rb,cmap,lo,hi):
    cm=colormaps[cmap]
    for c in chroms:
        for (st,v) in track[c]:
            s=th(c,st); e=th(c,min(st+WIN,chrlen[c]))
            ax.bar((s+e)/2,rb[1]-rb[0],width=abs(e-s),bottom=rb[0],color=cm(0.12+0.8*gnorm(v,lo,hi)),
                   edgecolor='none',align='center',zorder=2)
def ring_line(track,rb,color,lo,hi):
    for c in chroms:
        if not track[c]: continue
        ax.plot([th(c,st+WIN/2) for st,_ in track[c]],
                [rb[0]+(rb[1]-rb[0])*gnorm(v,lo,hi) for _,v in track[c]],color=color,lw=0.8,zorder=3)
ring_heat(gene,R_GENE,'Blues',g_lo,g_hi); ring_heat(te,R_TE,'OrRd',t_lo,t_hi); ring_line(gc,R_GC,C['green'],gc_lo,gc_hi)
chrom_cols=colormaps['twilight'](np.linspace(0.08,0.92,12))
for i,c in enumerate(chroms):
    s,e=spans[c]; mid=(s+e)/2
    ax.bar(mid,R_IDEO[1]-R_IDEO[0],width=abs(e-s),bottom=R_IDEO[0],color=chrom_cols[i],edgecolor='white',linewidth=0.6,zorder=4)
    deg=math.degrees(mid)%360; rot=deg-90 if deg<=180 else deg+90
    ax.text(mid,1.005,c.replace('Chr',''),rotation=rot,rotation_mode='anchor',ha='center',va='center',fontsize=7.2,fontweight='bold',color=C['ink'])
    # Mb ticks + labels every 10 Mb
    for mb in range(0,chrlen[c],10_000_000):
        t=th(c,mb); ax.plot([t,t],[R_IDEO[1],R_IDEO[1]+0.010],color='#444',lw=0.5,zorder=5)
        if mb>0: ax.text(t,R_IDEO[1]+0.026,f"{mb//1_000_000}",ha='center',va='center',fontsize=3.7,color='#777',rotation=rot)
    for end in (0,chrlen[c]):
        ax.plot(th(c,end),R_IDEO[1]+0.018,marker='o',ms=2.2,mfc=C['vermilion'],mec='none',zorder=6)
    # centromere proxy (max-TE window)
    if c in cen:
        ax.plot(th(c,cen[c]),(R_IDEO[0]+R_IDEO[1])/2,marker='v',ms=3.2,mfc='#222',mec='white',mew=0.3,zorder=7)
# pigmentation loci leader lines + labels on outer ring
for c,pos,lab in LOCI:
    t=th(c,pos)
    ax.plot([t,t],[R_IDEO[1]+0.04,R_IDEO[1]+0.128],color=C['purple'],lw=0.9,zorder=7)
    ax.plot(t,R_IDEO[1]+0.04,marker='*',ms=6,mfc=C['purple'],mec='white',mew=0.3,zorder=8)
    deg=math.degrees(t)%360; ha='left' if deg<180 else 'right'
    rot=deg-90 if deg<=180 else deg+90
    ax.text(t,R_IDEO[1]+0.142,lab,rotation=rot,rotation_mode='anchor',ha='center',va='center',fontsize=5.3,color=C['purple'],fontweight='bold')
# centre
ax.text(0.5,0.5,"ZN65\npurple rice (indica)\n\n395.1 Mb · 12 gap-free chr\n24 telomeres · QV 53.6\nLAI 15.5 · BUSCO 99.6%\n42,090 genes · 56.6% TE",
        transform=ax.transAxes,ha='center',va='center',fontsize=7.0,color=C['ink'],linespacing=1.3,
        bbox=dict(boxstyle='round,pad=0.6',fc='#F6F8FA',ec='#D9D9D9',lw=0.7))
fig.text(0.5,0.052,"Rings (outer to inner): ideogram (Mb ticks; stars, pigmentation loci; triangles, approx. centromere [max-TE window]; orange dots, telomeres)  |  gene density  |  TE coverage  |  GC",ha='center',fontsize=5.8,color=C['ink'])
fig.text(0.5,0.034,"Heatmaps: gene/TE per 500-kb window, globally 2–98%% percentile-scaled (comparable across chromosomes); bedtools coverage.",ha='center',fontsize=5.6,color='#777')
# real-unit colourbars
cax1=fig.add_axes([0.15,0.115,0.17,0.013]); cax2=fig.add_axes([0.45,0.115,0.17,0.013]); cax3=fig.add_axes([0.75,0.115,0.16,0.013])
cb1=mpl.colorbar.ColorbarBase(cax1,cmap=colormaps['Blues'],norm=Normalize(g_lo,g_hi),orientation='horizontal')
cb1.set_label('Gene density (genes / 500 kb)',fontsize=5.8); cb1.ax.tick_params(labelsize=5)
cb2=mpl.colorbar.ColorbarBase(cax2,cmap=colormaps['OrRd'],norm=Normalize(t_lo*100,t_hi*100),orientation='horizontal')
cb2.set_label('TE coverage (%)',fontsize=5.8); cb2.ax.tick_params(labelsize=5)
cax3.plot([0,1],[0.5,0.5],color=C['green'],lw=1.2); cax3.set_xlim(0,1); cax3.set_ylim(0,1); cax3.axis('off')
cax3.set_title(f"GC content ring: {gc_lo*100:.0f}–{gc_hi*100:.0f}%% (mean 43.7%%)",fontsize=5.6,pad=1)
fig.suptitle("Telomere-to-telomere genome landscape of purple rice ZN65",fontsize=11,y=0.965,fontweight='bold')
plt.savefig("Fig1_genome_landscape.png",bbox_inches='tight'); plt.savefig("Fig1_genome_landscape.pdf",bbox_inches='tight')
print(f"saved Fig1; gene {g_lo:.0f}-{g_hi:.0f}/500kb, TE {t_lo*100:.0f}-{t_hi*100:.0f}%, GC {gc_lo*100:.0f}-{gc_hi*100:.0f}%, cen markers {len(cen)}")
