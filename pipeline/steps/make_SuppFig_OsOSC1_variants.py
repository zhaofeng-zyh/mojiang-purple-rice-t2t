#!/usr/bin/env python3
"""Supplementary Fig. OsOSC1 (Os02g0139500 / ZN652G0336) variant map: ZN65 vs Nipponbare.
Reproducible from the MAFFT locus alignment (loci_pair.aln) + ZN65 exon coords (exons.txt).
Usage: python make_SuppFig_OsOSC1_variants.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pubstyle as nature_style
font, C = nature_style.apply()
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
from matplotlib.lines import Line2D

D = "/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/05_OsOSC1变异分析/"
# ZN65 locus window Chr2:2,059,636-2,072,871 (revcomp -> gene 5'->3'); gene TSS~2001, TES~11236
def rev(c): return 2072871 - c + 1
exons = []
for ln in open(D + "exons.txt"):
    s, e = map(int, ln.split())
    exons.append((rev(e), rev(s)))
exons.sort()
gs = min(s for s, e in exons); ge = max(e for s, e in exons)
Lwin = 13236
# large InDels (rev-locus coords) + type
indels = [(1, 1359, 'ZN65 insertion (MULE-MuDR TE)', C['red']),
          (10327, 52, 'ZN65 insertion (intron)', C['orange']),
          (13236, 132, 'Nipponbare insertion', C['grey'])]

fig, ax = plt.subplots(figsize=(8.0, 2.7))
ax.axvspan(1/1000, gs/1000, color='#FDF2E0', zorder=0)
ax.add_line(plt.Line2D([gs/1000, ge/1000], [0, 0], color='k', lw=1.0, zorder=1))
for s, e in exons:
    ax.add_patch(Rectangle((s/1000, -0.16), max((e-s)/1000, 0.03), 0.32, color=C['blue'], zorder=2))
ax.text((1+gs)/2/1000, -0.62, "5′ promoter", ha='center', fontsize=6.5, color='#9A6A00')
ax.text((gs+ge)/2/1000, -0.62, "OsOSC1 gene body (19 exons; CDS ~99.3% identical)", ha='center', fontsize=6.5, color=C['blue'])
def kb(n): return f"{n/1000:.2f} kb" if n>=1000 else f"{n} bp"
# promoter MULE-MuDR shown as a coloured BLOCK in the promoter (not an arrow)
pblk_e=gs; pblk_s=max(1,gs-1359)
ax.add_patch(Rectangle((pblk_s/1000,-0.17),(pblk_e-pblk_s)/1000,0.34,facecolor=C['red'],edgecolor='white',lw=0.4,zorder=4))
ax.text((pblk_s+pblk_e)/2/1000,0.40,"1.36-kb MULE-MuDR\ninsertion (ZN65)",ha='center',va='bottom',fontsize=6.0,color=C['red'],fontweight='bold')
# smaller InDels as size-labelled arrows
for pos, sz, lab, col in indels[1:]:
    ax.add_patch(FancyArrow(pos/1000, 0.22, 0, 0.18, width=0, head_width=0.16, head_length=0.1, color=col, zorder=4))
    ax.text(pos/1000, 0.44, f"{kb(sz)}\n{lab.split('(')[0].strip()}", ha='center', va='bottom', fontsize=5.4, color=col)
ax.legend(handles=[Rectangle((0,0),1,1,fc=C['red'],label='promoter TE block (MULE-MuDR)'),
                   Line2D([0],[0],marker='v',linestyle='None',markersize=6,markerfacecolor=C['orange'],markeredgecolor='none',label='ZN65 intron insertion'),
                   Line2D([0],[0],marker='v',linestyle='None',markersize=6,markerfacecolor=C['grey'],markeredgecolor='none',label='Nipponbare insertion')],
          loc='upper center', bbox_to_anchor=(0.5, -0.30), ncol=3, fontsize=5.8, frameon=False)
ax.set_xlim(-0.3, (Lwin+300)/1000); ax.set_ylim(-0.95, 0.95); ax.set_yticks([])
ax.set_xlabel("ZN65 OsOSC1 locus (Chr2:2,059,636–2,072,871; gene shown 5′→3′)  —  kb", fontsize=7.5)
ax.tick_params(labelsize=6.5)
for sp in ['top', 'right', 'left']: ax.spines[sp].set_visible(False)
ax.set_title("OsOSC1 (Os02g0139500) ZN65 vs Nipponbare: locus 98.93% identity, 64 SNPs, main SV = 1.36 kb promoter TE",
             fontsize=7.6)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "SuppFig_OsOSC1_variants")
plt.savefig(out + ".png"); plt.savefig(out + ".pdf")
print("saved", out + ".png/.pdf")
