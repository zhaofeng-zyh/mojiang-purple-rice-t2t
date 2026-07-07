#!/usr/bin/env python3
"""Supplementary Fig. Phylogenetic placement of the proposal gene Os02g0139500/OsOSC1
(a, OSC family) and the Chr2 MYB candidate ZN652G3275 (b, R2R3-MYB families).
Reproducible from IQ-TREE .treefile outputs. Usage: python make_SuppFig_phylogenies.py
"""
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pubstyle as nature_style
font, C = nature_style.apply()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def parse_newick(nwk):
    pos = [0]
    def parse():
        node = {'children': [], 'name': None, 'support': None, 'bl': 0.0}
        if nwk[pos[0]] == '(':
            pos[0] += 1
            while True:
                node['children'].append(parse())
                if nwk[pos[0]] == ',': pos[0] += 1; continue
                if nwk[pos[0]] == ')': pos[0] += 1; break
            m = re.match(r'([0-9.]+)?', nwk[pos[0]:])
            if m and m.group(1): node['support'] = m.group(1); pos[0] += m.end()
        m = re.match(r'([^,():;]+)?(?::(-?[0-9.eE]+))?', nwk[pos[0]:])
        if m:
            if m.group(1) and not node['children']: node['name'] = m.group(1)
            if m.group(2): node['bl'] = float(m.group(2))
            pos[0] += m.end()
        return node
    return parse()

def layout(root):
    leaves = []
    def collect(n):
        if not n['children']: leaves.append(n)
        else: [collect(c) for c in n['children']]
    collect(root)
    yo = {id(l): i for i, l in enumerate(leaves)}
    def xp(n, x=0.0):
        n['x'] = x + n['bl']
        if n['children']:
            for c in n['children']: xp(c, n['x'])
            n['y'] = sum(c['y'] for c in n['children']) / len(n['children'])
        else: n['y'] = yo[id(n)]
    xp(root)
    return leaves

def draw(ax, root, colmap, pretty, bold_set, xmax):
    def rec(n):
        if n['children']:
            ys = [c['y'] for c in n['children']]
            ax.plot([n['x'], n['x']], [min(ys), max(ys)], '-', color='k', lw=0.7)
            for c in n['children']:
                nm = c.get('name')
                col = colmap(nm) if not c['children'] else 'k'
                ax.plot([n['x'], c['x']], [c['y'], c['y']], '-', color=col,
                        lw=1.6 if nm in bold_set else 0.8)
                rec(c)
            if n['support']:
                ax.text(n['x'] - xmax*0.02, n['y'] + 0.42, n['support'], fontsize=5.0, ha='right', va='bottom', color='#666')
        else:
            nm = n['name']
            ax.text(n['x'] + xmax*0.025, n['y'], pretty.get(nm, nm), va='center',
                    fontsize=6.0, color=colmap(nm), fontweight=('bold' if nm in bold_set else 'normal'))
    rec(root)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.4,5.6))

# ---- panel a: OSC family ----
oscf = "/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/04_系统发育/OSC_tree.treefile"
rootA = parse_newick(open(oscf).read().strip()); leavesA = layout(rootA)
hotA = {'ZN652G0336_Chr2', 'Os02g0139500_OsOSC1_Nipponbare'}
prettyA = {'ZN652G0336_Chr2':'ZN652G0336 (ZN65)  ★','Os02g0139500_OsOSC1_Nipponbare':'Os02g0139500 / OsOSC1  ★',
 'ZN652G0345_Chr2_Os02g0139500_ortholog':'ZN652G0345 (ZN65)','Os02g0139700_CAS_Nipponbare':'Os02g0139700 / CAS',
 'CAS_rice_cycloartenol_Q6Z2X6':'Cycloartenol synthase, rice','CAS1_Arabidopsis_cycloartenol_P38605':'Cycloartenol synthase, Arabidopsis',
 'AchilleolB_synthase_rice_Q2R712':'Achilleol B synthase, rice','Parkeol_synthase_rice_H2KWF1':'Parkeol synthase, rice',
 'Lanosterol_synthase_human_OUTGROUP_P48449':'Lanosterol synthase, human (outgroup)',
 'ZN652G0349_Chr2':'ZN652G0349 (ZN65)','ZN652G0351_Chr2':'ZN652G0351 (ZN65)','ZN6511G1908_Chr11':'ZN6511G1908 (ZN65)',
 'ZN6511G1142_Chr11':'ZN6511G1142 (ZN65)','ZN6511G0633_Chr11':'ZN6511G0633 (ZN65)'}
xmaxA = max(n['x'] for n in leavesA)
draw(axA, rootA, lambda nm: C['red'] if nm in hotA else 'k', prettyA, hotA, xmaxA)
axA.set_xlim(-0.02, xmaxA*1.9); axA.set_ylim(-0.6, len(leavesA)-0.4); axA.set_yticks([])
axA.set_xlabel("substitutions / site", fontsize=7); axA.tick_params(labelsize=6)
for s in ['top','right','left']: axA.spines[s].set_visible(False)
axA.set_title("a  Os02g0139500 / OsOSC1 within oxidosqualene cyclases", fontsize=7.6, loc='left', fontweight='bold')

# ---- panel b: MYB ----
mybf = "/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/06_MYB_ZN652G3275/MYB_tree.treefile"
rootB = parse_newick(open(mybf).read().strip()); leavesB = layout(rootB)
antho = {'OsC1_Os06g0205100_anthocyanin_MYB','ZmC1_maize_anthocyanin','AtPAP1_MYB75_anthocyanin','AtPAP2_MYB90_anthocyanin','AtMYB113_anthocyanin'}
lignin = {'ZN652G3275_ZN65_Chr2_QUERY','OsMYB58_63L_Os02g0695200','AtMYB58_lignin_SCW','AtMYB63_lignin_SCW'}
prettyB = {'ZN652G3275_ZN65_Chr2_QUERY':'ZN652G3275 (ZN65)  ★','OsMYB58_63L_Os02g0695200':'Os02g0695200 / OsMYB58-63',
 'AtMYB58_lignin_SCW':'AtMYB58 (lignin)','AtMYB63_lignin_SCW':'AtMYB63 (lignin)','AtMYB46_SCW_master':'AtMYB46 (SCW)',
 'AtMYB12_flavonol':'AtMYB12 (flavonol)','OsC1_Os06g0205100_anthocyanin_MYB':'OsC1 (anthocyanin)','ZmC1_maize_anthocyanin':'ZmC1 (anthocyanin)',
 'AtPAP1_MYB75_anthocyanin':'AtPAP1 (anthocyanin)','AtPAP2_MYB90_anthocyanin':'AtPAP2 (anthocyanin)','AtMYB113_anthocyanin':'AtMYB113 (anthocyanin)'}
def colB(nm): return C['red'] if nm in lignin else (C['blue'] if nm in antho else C['grey'])
xmaxB = max(n['x'] for n in leavesB)
draw(axB, rootB, colB, prettyB, lignin, xmaxB)
axB.set_xlim(-0.05, xmaxB*1.75); axB.set_ylim(-0.6, len(leavesB)-0.4); axB.set_yticks([])
axB.set_xlabel("substitutions / site", fontsize=7); axB.tick_params(labelsize=6)
for s in ['top','right','left']: axB.spines[s].set_visible(False)
axB.set_title("b  ZN652G3275 clusters with lignin MYBs, not anthocyanin MYBs", fontsize=7.6, loc='left', fontweight='bold')
axB.legend(handles=[Line2D([0],[0],color=C['red'],lw=2,label='lignin MYB (incl. ZN652G3275)'),
                    Line2D([0],[0],color=C['blue'],lw=2,label='anthocyanin MYB')],
           loc='lower right', fontsize=5.6, frameon=False)

fig.suptitle("Supplementary Figure: the proposal gene 'Os02g' is OsOSC1 (a triterpene cyclase), not an anthocyanin gene",
             fontsize=8, y=1.0)
plt.tight_layout(rect=[0,0,1,0.96])
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "SuppFig_Os02g_phylogenies")
plt.savefig(out + ".png"); plt.savefig(out + ".pdf")
print("saved", out + ".png/.pdf")
