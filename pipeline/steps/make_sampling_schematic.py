import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['font.sans-serif']=['Noto Sans CJK JP','Droid Sans Fallback','DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus']=False
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Ellipse, FancyArrow
fig=plt.figure(figsize=(9.2,8.0))
gs=fig.add_gridspec(2,1,height_ratios=[0.95,1.05],hspace=0.32)

# ===== Panel A: 籽粒结构与取样部位 =====
axA=fig.add_subplot(gs[0]); axA.set_xlim(0,10); axA.set_ylim(0,5); axA.axis('off')
axA.text(0.1,4.7,"A  水稻籽粒结构与取样部位",fontsize=12,fontweight='bold')
gx,gy=2.2,2.3
axA.add_patch(Ellipse((gx,gy),3.0,1.7,facecolor='#7A3B8F',edgecolor='#4A2356',lw=1.2,zorder=2))
axA.add_patch(Ellipse((gx,gy),2.78,1.5,facecolor='#C9A24B',edgecolor='none',zorder=3))
axA.add_patch(Ellipse((gx,gy),2.55,1.32,facecolor='#F3EFE2',edgecolor='none',zorder=4))
axA.add_patch(Ellipse((gx-1.05,gy-0.18),0.55,0.62,facecolor='#E2C275',edgecolor='#9A7A2A',lw=0.8,zorder=5))
axA.text(gx-1.05,gy-0.18,"胚",ha='center',va='center',fontsize=8,zorder=6)
axA.text(gx,gy,"胚乳\n(淀粉)",ha='center',va='center',fontsize=8,color='#777',zorder=6)
axA.annotate("颖壳/谷壳\n(已脱去)",xy=(gx+1.4,gy+0.7),xytext=(gx-0.3,4.4),fontsize=7.5,color='#999',ha='center',
             arrowprops=dict(arrowstyle='-',color='#bbb',lw=0.6,ls='--'))
axA.text(gx,0.45,"糙米纵切（颖壳已脱）",ha='center',fontsize=8,color='#555')
zx=6.0; bars=[("果皮 pericarp",'#7A3B8F',0.34),("种皮 testa",'#9B6FB0',0.16),("糊粉层 aleurone",'#C9A24B',0.30),("胚乳 endosperm",'#F3EFE2',1.5)]
y=4.0
for lab,col,w in bars:
    axA.add_patch(Rectangle((zx,y-w),2.3,w,facecolor=col,edgecolor='white',lw=0.6,zorder=3))
    axA.text(zx+2.45,y-w/2,lab,va='center',fontsize=8.2)
    y-=w
top=4.0; bot=4.0-(0.34+0.16+0.30)
axA.annotate('',xy=(zx-0.18,top),xytext=(zx-0.18,bot),arrowprops=dict(arrowstyle='|-|',color='#B0271A',lw=1.4))
axA.text(zx-0.35,(top+bot)/2,"米糠层 (bran)\n果皮+种皮+糊粉层\n花青素富集 → 取样部位",ha='right',va='center',fontsize=8.0,color='#B0271A',fontweight='bold')
axA.text(zx+1.15,4.35,"外层放大（横切层栈）",ha='center',fontsize=8,color='#555')
axA.annotate('',xy=(zx-0.05,3.2),xytext=(gx+1.45,gy+0.2),arrowprops=dict(arrowstyle='-|>',color='#999',lw=0.9,connectionstyle="arc3,rad=-0.15"))

# ===== Panel B: 发育期取样时间轴 + 中文分期 =====
axB=fig.add_subplot(gs[1]); axB.set_xlim(-3.5,34); axB.set_ylim(0,5); axB.axis('off')
axB.text(-3.5,4.75,"B  果皮发育期取样时间轴（以开花/抽穗为 0 DAF；中文分期）",fontsize=12,fontweight='bold')
TY=2.15   # timeline y
axB.add_line(plt.Line2D([0,33],[TY,TY],color='#333',lw=1.2))
stages=[(0,'#F6F3EA','0\n开花'),(5,'#F3E6D8','5'),(10,'#E9C6C0','10'),(15,'#C98FB0','15'),(20,'#8F4E8F','20'),(25,'#6A2F75','25'),(30,'#4A2356','30')]
for d,col,lab in stages:
    axB.add_patch(Ellipse((d,TY),1.7,0.82,facecolor=col,edgecolor='#555',lw=0.6,zorder=3))
    axB.add_line(plt.Line2D([d,d],[TY,TY-0.45],color='#333',lw=0.6))
    axB.text(d,TY-0.55,lab,ha='center',va='top',fontsize=7.6)
# 中文分期带
phases=[(3,13,'乳熟期','#FBEAD2','3–13 DAF'),(13,23,'蜡熟期','#E9C9DC','13–23 DAF'),(23,33,'完熟期','#CDB3DC','23–35 DAF')]
by0,by1=0.55,1.10
for a,b,name,col,rng in phases:
    axB.add_patch(FancyBboxPatch((a,by0),b-a-0.15,by1-by0,boxstyle="round,pad=0.02",
                  facecolor=col,edgecolor='#888',lw=0.6,zorder=2))
    axB.text((a+b)/2,(by0+by1)/2+0.12,name,ha='center',va='center',fontsize=9,fontweight='bold',color='#5A2A6A')
    axB.text((a+b)/2,(by0+by1)/2-0.16,rng,ha='center',va='center',fontsize=7,color='#7A5A88')
axB.text(16.5,0.18,"开花后天数 DAF（籽粒/果皮颜色随发育加深）",ha='center',fontsize=8.5,color='#555')
# 花青素起始
axB.annotate("花青素约 10–15 DAF 起在果皮积累",xy=(12,TY+0.55),xytext=(11.5,4.0),fontsize=8,color='#7A3B8F',ha='center',
             arrowprops=dict(arrowstyle='-|>',color='#7A3B8F',lw=1.0))
# G2 取样箭头
for d in [5,10,15,20]:
    axB.add_patch(FancyArrow(d,3.30,0,-0.55,width=0.0,head_width=0.7,head_length=0.25,color='#1A7A4A',zorder=5))
axB.text(12.5,3.52,"G2 果皮取样：5 / 10 / 15 / 20 DAF（各 ≥3 重复）+ 营养期对照",ha='center',fontsize=8.4,color='#1A7A4A',fontweight='bold')
# G1
axB.annotate("G1 花青素定量：完熟期糙米 → 碾出米糠层",xy=(30,1.15),xytext=(26.5,4.15),fontsize=8.2,color='#B0271A',ha='center',fontweight='bold',
             arrowprops=dict(arrowstyle='-|>',color='#B0271A',lw=1.0))
fig.suptitle("ZN65 果皮花青素(G1)/发育期表达(G2) —— 取样部位与时间示意",fontsize=11,y=0.975,fontweight='bold')
out="/mnt/data2/墨江紫米研究/13_Paper2_机制研究方案/SOP_取样示意图.png"
plt.savefig(out,dpi=200,bbox_inches='tight'); plt.savefig(out.replace('.png','.pdf'),bbox_inches='tight')
print("saved",out)
