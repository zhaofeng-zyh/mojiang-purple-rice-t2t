import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
from collections import defaultdict
W="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/11_基因组景观/"
def load(f,col):
    d=defaultdict(list)
    for ln in open(W+f):
        p=ln.split("\t")
        try: d[p[0]].append((int(p[1]),float(p[col])))
        except: pass
    return d
gene=load("gene_density.txt",3)
te=load("te_cov.txt",6)   # fraction covered (last col of bedtools coverage, 3-col -a -> col index 6)
gc=load("gc.txt",3)
chroms=[f"Chr{i}" for i in range(1,13)]
fig,axes=plt.subplots(12,1,figsize=(11,9),sharex=True)
for ax,ch in zip(axes,chroms):
    g=sorted(gene.get(ch,[])); t=sorted(te.get(ch,[]))
    if g:
        x=[p[0]/1e6 for p in g]; y=[p[1] for p in g]
        ax.fill_between(x,0,y,color="#1f6fb4",alpha=0.75,lw=0,label="gene density")
        ymax=max(y) if y else 1
        ax.set_ylim(0,ymax*1.05)
    if t:
        xt=[p[0]/1e6 for p in t]; yt=[p[1]*ax.get_ylim()[1] for p in t]
        ax.fill_between(xt,0,yt,color="#c0392b",alpha=0.45,lw=0,label="TE fraction")
    ax.set_ylabel(ch,rotation=0,ha="right",va="center",fontsize=8)
    ax.set_yticks([]); 
    for s in ["top","right","left"]: ax.spines[s].set_visible(False)
axes[0].legend(loc="upper right",fontsize=8,ncol=2,frameon=False)
axes[-1].set_xlabel("Chromosome position (Mb)")
fig.suptitle("Genome landscape of purple rice ZN65 (T2T): gene density (blue) vs transposon density (red)\n"
             "12 gap-free chromosomes, 395 Mb, 42,090 genes, 56.6% repeats — genes on arms, TEs pericentromeric",fontsize=11)
plt.tight_layout(rect=[0,0,1,0.95]); plt.savefig(W+"genome_landscape.png",dpi=170)
print("saved genome_landscape.png")
