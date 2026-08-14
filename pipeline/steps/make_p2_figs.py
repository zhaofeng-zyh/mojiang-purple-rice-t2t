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
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams['font.sans-serif']=['Noto Sans CJK JP','Droid Sans Fallback','DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus']=False
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow, FancyBboxPatch, Polygon
D="/mnt/data2/墨江紫米研究/13_Paper2_机制研究方案/"

# ============ FIG A: OsB2 位点 + 分子试剂图 ============
fig=plt.figure(figsize=(10.2,7.4)); gs=fig.add_gridspec(2,1,height_ratios=[1.25,1.0],hspace=0.42)
ax=fig.add_subplot(gs[0])
x0,x1=27995.5,28029.0   # kb
ax.set_xlim(x0,x1); ax.set_ylim(-3.2,3.2); ax.axis('off')
def kb(b): return b/1000.0
ax.text(x0,3.0,"A  OsB2/Kala4 位点（Chr4，− 链）与分子试剂设计",fontsize=12,fontweight='bold')
# baseline
ax.add_line(plt.Line2D([x0,x1],[0,0],color='#999',lw=0.8,zorder=1))
# gene body (mRNA 27996908-28021147, -strand → 5'在右)
ax.add_line(plt.Line2D([kb(27996908),kb(28021147)],[0,0],color='#333',lw=1.4,zorder=2))
exons=[(28020506,28020642),(28013739,28013993),(28012980,28013076),(28012729,28012743),(28012559,28012615),(27997503,27998159),(27997288,27997425)]
for s,e in exons:
    ax.add_patch(Rectangle((kb(s),-0.28),kb(e)-kb(s),0.56,facecolor='#2E5A87',edgecolor='k',lw=0.4,zorder=4))
ax.annotate('',xy=(kb(27996908),0),xytext=(kb(27998000),0),arrowprops=dict(arrowstyle='-|>',color='#333',lw=1.2))  # 5'->3' left
ax.text(kb(28009000),0.55,"OsB2 ZN654G2687 (7 exons, bHLH)",ha='center',fontsize=8,color='#2E5A87',fontweight='bold')
# intronic TEs
ax.add_patch(Rectangle((kb(27999746),-0.18),kb(28012241)-kb(27999746),0.36,facecolor='#C45B5B',alpha=0.45,edgecolor='none',zorder=3))
ax.text(kb(28006000),-0.55,"Gypsy RETRO2B 12.5kb (intron)",ha='center',fontsize=7,color='#A03030')
ax.add_patch(Rectangle((kb(28014081),-0.18),kb(28020405)-kb(28014081),0.36,facecolor='#E0A458',alpha=0.5,edgecolor='none',zorder=3))
ax.text(kb(28017200),-0.55,"LINE1 6.3kb (intron)",ha='center',fontsize=7,color='#9A6A20')
# promoter region
ax.add_line(plt.Line2D([kb(28021147),kb(28028800)],[0,0],color='#7AA0C0',lw=1.0,ls=':',zorder=2))
ax.text(kb(28021300),0.55,"启动子(5'上游)",ha='left',fontsize=7.5,color='#5A80A0')
# ZN65-specific insertion 28022131-28028000
i0,i1=28022131,28028000
ax.add_patch(FancyBboxPatch((kb(i0),0.9),kb(i1)-kb(i0),0.7,boxstyle="round,pad=0.02",facecolor='#B0271A',edgecolor='#6A1810',lw=0.8,zorder=5))
ax.text(kb((i0+i1)/2),1.25,f"ZN65 特异启动子插入 ~5.9 kb\n(SINE/Helitron/hAT/PIF 簇)",ha='center',va='center',fontsize=7.4,color='white',fontweight='bold')
ax.annotate('',xy=(kb((i0+i1)/2),0.05),xytext=(kb((i0+i1)/2),0.88),arrowprops=dict(arrowstyle='-',color='#B0271A',lw=0.8,ls='--'))
# sgRNA at exon2 (28013739-993)
sg=kb((28013739+28013993)/2)
ax.annotate("sgRNA×2\n(exon2, 唯一)",xy=(sg,0.30),xytext=(sg,2.35),ha='center',fontsize=7.2,color='#1A7A4A',fontweight='bold',
            arrowprops=dict(arrowstyle='-|>',color='#1A7A4A',lw=1.1))
# genotyping primers (positions approx)
prim=[("F-flank5",28021650,'+','#333',-1.62,'right'),("R-ins5",28022600,'-','#B0271A',-2.02,'left'),
      ("F-ins3",28027450,'+','#B0271A',-1.62,'right'),("R-flank3",28028450,'-','#333',-2.02,'left')]
for nm,pos,strand,col,ly,ha in prim:
    dx=0.16 if strand=='+' else -0.16
    ax.add_patch(FancyArrow(kb(pos)-dx/2,-1.30,dx,0,width=0.0,head_width=0.16,head_length=0.11,color=col,zorder=6))
    ax.add_line(plt.Line2D([kb(pos),kb(pos)],[-1.30,ly+0.10],color=col,lw=0.4,ls=':'))
    ax.text(kb(pos),ly,nm,ha=ha,va='top',fontsize=6.4,color=col)
ax.text(kb(28025000),-2.62,"断点分型引物（共显性 存在/缺失）",ha='center',fontsize=7.2,color='#555')
# Mb ticks
for t in range(27996,28029,4):
    ax.add_line(plt.Line2D([kb(t*1000)]*2,[-0.06,0.06],color='#999',lw=0.6))
    ax.text(kb(t*1000),-0.95 if t<27999 else -3.0,f"{t/1000:.3f}",ha='center',fontsize=5.6,color='#888')
ax.text(x1,-3.0,"Mb",fontsize=6,color='#888')

# ---- panel B: 分型胶图判读 ----
axb=fig.add_subplot(gs[1]); axb.set_xlim(0,10); axb.set_ylim(0,6); axb.axis('off')
axb.text(0,5.7,"B  插入断点 PCR 分型判读（共分离/种质 panel）",fontsize=11,fontweight='bold')
axb.text(0.2,5.15,"反应①: F-flank5 + R-ins5  →  仅\"含插入\"等位出带（紫连接5')",fontsize=8,color='#B0271A')
axb.text(0.2,4.78,"反应②: F-ins3 + R-flank3   →  仅\"含插入\"等位出带（紫连接3')",fontsize=8,color='#B0271A')
axb.text(0.2,4.41,"反应③: F-flank5 + R-flank3 →  仅\"无插入\"等位出小带（白；含插入太长~5.9kb不扩）",fontsize=8,color='#333')
# 三个泳道示意
lanes=[("紫 纯合\n(ins/ins)",[1,1,0]),("杂合\n(ins/–)",[1,1,1]),("白 纯合\n(–/–)",[0,0,1])]
lx=[2.0,5.0,8.0]
for (lab,bands),cx in zip(lanes,lx):
    axb.add_patch(Rectangle((cx-0.85,0.35),1.7,2.6,facecolor='#F3F3F0',edgecolor='#aaa',lw=0.6))
    axb.text(cx,3.18,lab,ha='center',fontsize=8,fontweight='bold')
    ys=[2.55,2.0,0.95]; labs=['①','②','③']; cols=['#B0271A','#B0271A','#333']
    for y,b,l,c in zip(ys,bands,labs,cols):
        axb.text(cx-1.05,y,l,ha='right',fontsize=7,color=c)
        if b: axb.add_patch(Rectangle((cx-0.6,y-0.06),1.2,0.16,facecolor=c,edgecolor='none'))
        else: axb.text(cx,y,"—",ha='center',va='center',fontsize=8,color='#bbb')
fig.suptitle("Paper2 图：OsB2 位点分子试剂与基因型分型",fontsize=11,y=0.98,fontweight='bold')
plt.savefig(D+"Paper2_Fig_OsB2位点与分子试剂.png",dpi=200,bbox_inches='tight')
plt.savefig(D+"Paper2_Fig_OsB2位点与分子试剂.pdf",bbox_inches='tight'); plt.close()
print("saved Fig A (locus)")

# ============ FIG C: 功能验证证据链与实验策略 ============
fig2,ax=plt.subplots(figsize=(10.2,5.6)); ax.set_xlim(0,10); ax.set_ylim(0,6); ax.axis('off')
ax.text(0,5.7,"功能验证策略：基因型 → 表达 → 代谢物 → 表型（因果闭环）",fontsize=12,fontweight='bold')
boxes=[("基因型\n启动子TE插入/SV",0.3,'#D9C2E0'),("转录\nOsB2/MBW/结构基因",2.7,'#C2D9E0'),("代谢物\n花青素 C3G/P3G",5.1,'#E0D2C2'),("表型\n果皮紫色",7.5,'#E0C2C8')]
for lab,x,col in boxes:
    ax.add_patch(FancyBboxPatch((x,3.7),2.1,1.0,boxstyle="round,pad=0.04",facecolor=col,edgecolor='#666',lw=0.8))
    ax.text(x+1.05,4.2,lab,ha='center',va='center',fontsize=8.5,fontweight='bold')
for x in [2.55,4.95,7.35]:
    ax.annotate('',xy=(x+0.18,4.2),xytext=(x-0.05,4.2),arrowprops=dict(arrowstyle='-|>',color='#444',lw=1.4))
# assays under each linkage
assays=[("Aim3 编辑/启动子互换\nCRISPR-KO·过表达·promoter-swap",1.35,'#7A3B8F'),
        ("Aim2 转录组+qPCR\n(比对ZN65自有基因组)",3.75,'#2E6A8F'),
        ("Aim2 靶向代谢组\nLC-MS/MS MRM 绝对定量",6.15,'#9A6A20'),
        ("Aim1 表型/显微定位\nHPLC·切片·LCM",8.55,'#A03048')]
for lab,x,col in assays:
    ax.annotate('',xy=(x,3.65),xytext=(x,3.05),arrowprops=dict(arrowstyle='-',color=col,lw=0.8,ls=':'))
    ax.add_patch(FancyBboxPatch((x-1.0,1.9),2.0,1.05,boxstyle="round,pad=0.03",facecolor='white',edgecolor=col,lw=1.0))
    ax.text(x,2.42,lab,ha='center',va='center',fontsize=7.0,color=col)
# mechanism layer
ax.add_patch(FancyBboxPatch((1.2,0.4),7.6,1.0,boxstyle="round,pad=0.04",facecolor='#F3EEF6',edgecolor='#7A3B8F',lw=1.0))
ax.text(5.0,0.9,"机制层(Aim3)：MBW 复合体  OsB2(bHLH)×OsC1(MYB)×OsTTG1(WD40)  → 反式激活 DFR/ANS/3GT 启动子\n"
        "Y2H·BiFC·Split-LUC(互作) ｜ EMSA·Y1H·ChIP-qPCR·dual-LUC(结合/激活)",
        ha='center',va='center',fontsize=7.6,color='#5A2A6A')
ax.annotate('',xy=(5.0,1.45),xytext=(5.0,1.85),arrowprops=dict(arrowstyle='-|>',color='#7A3B8F',lw=1.2))
plt.savefig(D+"Paper2_Fig_功能验证策略.png",dpi=200,bbox_inches='tight')
plt.savefig(D+"Paper2_Fig_功能验证策略.pdf",bbox_inches='tight'); plt.close()
print("saved Fig C (strategy)")
