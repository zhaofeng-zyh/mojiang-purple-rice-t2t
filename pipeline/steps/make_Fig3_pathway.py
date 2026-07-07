import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import pubstyle; font,C=pubstyle.apply()
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch

TSV="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/09_花青素通路清单/ZN65_花青素通路基因清单.tsv"
RBH="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/18_通路拷贝数_A11/best.tsv"
abbr={"PAL":"PAL","C4H":"C4H","4CL":"4CL","CHS":"CHS","CHI":"CHI","F3H":"F3H","F3'H":"F3'H",
      "F3'5'H":"F3'5'H","DFR":"DFR","ANS/LDOX":"ANS","UFGT/3GT":"3GT","LAR":"LAR","ANR":"ANR","FLS":"FLS"}
zn2nip={}
try:
    for ln in open(RBH):
        f=ln.rstrip("\n").split("\t")
        if len(f)>=2: zn2nip[f[0].split(".")[0]]=f[1]
except FileNotFoundError: pass
cnt={v:0 for v in abbr.values()}; nipset={v:set() for v in abbr.values()}
for ln in open(TSV):
    f=ln.rstrip("\n").split("\t"); step=f[0].strip(); key=step.split()[0] if step else ""
    if key in abbr:
        a=abbr[key]; cnt[a]+=1
        if len(f)>2 and f[2].strip() in zn2nip: nipset[a].add(zn2nip[f[2].strip()])
nip={a:len(s) for a,s in nipset.items()}
def cn(a): return f"{cnt.get(a,0)} / {nip.get(a,0)}"
fh=cnt["F3'H"]; fhN=nip["F3'H"]; ffh=cnt["F3'5'H"]; ffhN=nip["F3'5'H"]

fig,ax=plt.subplots(figsize=(9.0,5.8)); ax.set_xlim(-0.4,13.7); ax.set_ylim(-3.2,3.0); ax.axis('off')
def box(x,y,label,a,kind='enz',w=0.96,h=0.7):
    col={'enz':'#E8EEF5','commit':'#FAD9D2','term':'#E7D6F0','reg':'#DFF2EA'}[kind]
    ec={'enz':C['blue'],'commit':C['vermilion'],'term':C['purple'],'reg':C['green']}[kind]
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.03,rounding_size=0.10',fc=col,ec=ec,lw=1.1,zorder=3))
    ax.text(x,y+0.11,label,ha='center',va='center',fontsize=8.2,fontweight='bold',color=C['ink'],zorder=4)
    if a is not None: ax.text(x,y-0.17,cn(a),ha='center',va='center',fontsize=6.2,color=ec,zorder=4)
def arrow(x1,y1,x2,y2,color='#666',lw=1.3,style='-|>',rad=0.0):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle=style,mutation_scale=9,color=color,lw=lw,zorder=2,
                 connectionstyle=f"arc3,rad={rad}"))
chain=[("PAL","PAL",'enz'),("C4H","C4H",'enz'),("4CL","4CL",'enz'),("CHS","CHS",'enz'),("CHI","CHI",'enz'),
       ("F3H","F3H",'enz'),("DFR","DFR",'commit'),("ANS","ANS",'commit'),("3GT","3GT",'commit')]
xs=[0.85+i*1.24 for i in range(len(chain))]; y0=1.15
for (lab,a,kind),x in zip(chain,xs): box(x,y0,lab,a,kind)
for i in range(len(xs)-1): arrow(xs[i]+0.49,y0,xs[i+1]-0.49,y0)
# metabolites — well below the chain, tiny, sparse
for x,txt in [((xs[3]+xs[4])/2,"naringenin\nchalcone"),((xs[5]+xs[6])/2,"DHK/DHQ"),((xs[6]+xs[7])/2,"leuco-\nanthocyanidin"),((xs[7]+xs[8])/2,"anthocyanidin")]:
    ax.text(x,y0+0.55,txt,ha='center',va='bottom',fontsize=5.0,color='#9aa',style='italic')
xt=xs[-1]+1.55; box(xt,y0,"Anthocyanins",None,'term',w=1.7,h=0.74)
ax.text(xt,y0-0.22,"C3G · P3G",ha='center',fontsize=6.0,color=C['purple'])
ax.text(xt,y0+0.62,"to be quantified\n(HPLC-DAD, G1)",ha='center',va='bottom',fontsize=5.6,color=C['red'],style='italic')
arrow(xs[-1]+0.49,y0,xt-0.86,y0)
ax.text(-0.3,y0,"Phe",ha='center',va='center',fontsize=7.2,color='#555'); arrow(-0.16,y0,xs[0]-0.49,y0)
# upper branches
box(xs[5]+0.2,2.55,"F3'H / F3'5'H",None,'enz',w=1.45,h=0.58); ax.text(xs[5]+0.2,2.34,f"{fh}/{fhN} . {ffh}/{ffhN}",ha="center",fontsize=5.6,color=C["blue"])
arrow(xs[5]+0.2,2.26,xs[5]+0.15,y0+0.37,color=C['blue'],lw=1.0)
box(xs[4]-0.35,2.55,"FLS","FLS",'enz',w=0.86,h=0.58)
arrow(xs[4]-0.35,y0+0.37,xs[4]-0.35,2.26,color='#999',lw=1.0)
ax.text(xs[4]-0.95,2.55,"flavonols",ha='right',va='center',fontsize=5.8,color='#888')
# proanthocyanidin branch — LEFT-lower, off DFR, clear of MBW arrows
lx,ly=xs[5]-0.05,-1.45
box(lx,ly,"LAR / ANR",None,'enz',w=1.5,h=0.58); ax.text(lx,ly-0.20,f"{cnt['LAR']}/{nip['LAR']} · {cnt['ANR']}/{nip['ANR']}",ha='center',fontsize=5.6,color=C['blue'])
arrow(xs[6]-0.3,y0-0.37,lx+0.4,ly+0.30,color='#999',lw=1.0,rad=0.2)
ax.text(lx-0.82,ly,"proanthocyanidins\n(Rc branch)",ha='right',va='center',fontsize=5.8,color='#888')
# MBW — bottom, wide; VERTICAL non-crossing arrows to committed steps (right side, clear of LAR/ANR)
mleft,mright=xs[7]-0.55,xs[8]+0.55; mx=(mleft+mright)/2; my=-2.45
ax.add_patch(FancyBboxPatch((mleft-0.4,my-0.36),(mright-mleft)+0.8,0.72,boxstyle='round,pad=0.03,rounding_size=0.12',fc='#DFF2EA',ec=C['green'],lw=1.2,zorder=3))
ax.text(mx,my+0.13,"MBW complex  (model)",ha='center',fontsize=7.6,fontweight='bold',color=C['green'])
ax.text(mx,my-0.16,"OsC1 (MYB) · OsB2 (bHLH) · OsTTG1 (WD40)",ha='center',fontsize=6.0,color=C['ink'])
for x in (xs[7],xs[8]):           # vertical arrows to ANS, 3GT (DFR fed via leuco-branch already)
    arrow(x,my+0.37,x,y0-0.37,color=C['green'],lw=1.0,style='-|>')
arrow(xs[6],my+0.30,xs[6],y0-0.40,color=C['green'],lw=1.0,style='-|>',rad=0.0)
ax.text(mright+0.7,my,"transcriptional\nactivation (model)",ha='left',va='center',fontsize=5.8,color=C['green'])
# footnote
ax.text(-0.3,-3.05,"× = ZN65 / Nipponbare copies (KEGG KO, reciprocal-best-hit; may include paralogues/pseudogenes — curation pending).   "
        "SV-bearing loci: Kala4/OsB2, Kala1/OsDFR, Rc;  OsC1/Kala3 conserved.",fontsize=5.6,color='#666',ha='left')
ax.legend(handles=[Patch(fc='#E8EEF5',ec=C['blue'],label='enzyme (ZN65/NIP copies)'),
                   Patch(fc='#FAD9D2',ec=C['vermilion'],label='anthocyanin-committed'),
                   Patch(fc='#E7D6F0',ec=C['purple'],label='product (not yet quantified)'),
                   Patch(fc='#DFF2EA',ec=C['green'],label='MBW (model)')],
          loc='upper center',ncol=4,frameon=False,fontsize=6.2,bbox_to_anchor=(0.5,1.10))
ax.set_title("Flavonoid/anthocyanin pathway in ZN65 vs Nipponbare: complete, conserved enzyme complement; pigmentation set by MBW regulators",fontsize=8.4,y=1.15,fontweight='bold')
plt.savefig("Fig3_pathway.png",bbox_inches='tight',dpi=200); plt.savefig("Fig3_pathway.pdf",bbox_inches='tight')
print("saved Fig3_pathway v3")
