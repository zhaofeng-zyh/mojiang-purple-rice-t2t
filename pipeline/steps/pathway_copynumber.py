import collections
BASE="/mnt/data2/墨江紫米研究"
OGF="orthofinder_in/OrthoFinder/Results_Jul05/Orthogroups/Orthogroups.txt"
PW=f"{BASE}/07_分析_Os02g基因鉴定/09_花青素通路清单/ZN65_花青素通路基因清单.tsv"
gene2og={}; og_nipgenes=collections.defaultdict(set); og_zngenes=collections.defaultdict(set)
for line in open(OGF):
    if ':' not in line: continue
    og,rest=line.split(':',1)
    for g in rest.split():
        if g.startswith('ZN6'):
            gid=g.split('.')[0]; gene2og[gid]=og; og_zngenes[og].add(gid)
        elif g.startswith('Os'):
            og_nipgenes[og].add(g)
step_genes=collections.defaultdict(list); step_ko={}; order=[]
hdr=True
for line in open(PW):
    if hdr: hdr=False; continue
    f=line.rstrip('\n').split('\t')
    if len(f)<3: continue
    step,ko,gene=f[0],f[1],f[2].split('.')[0]
    if step not in order: order.append(step)
    step_genes[step].append(gene); step_ko[step]=ko
print(f"{'Step':<12}{'KO':<9}{'ZN65_KEGG':>10}{'ZN65_OG':>9}{'NIP_OG':>8}   (OG-level = distinct genes in spanned orthogroups)")
rows=[]; t1=t2=t3=0
for step in order:
    genes=set(step_genes[step])
    ogs=set(gene2og[g] for g in genes if g in gene2og)
    zn_kegg=len(genes)
    zn_og=len(set().union(*[og_zngenes[o] for o in ogs])) if ogs else 0
    nip_og=len(set().union(*[og_nipgenes[o] for o in ogs])) if ogs else 0
    en=step.split()[0]
    print(f"{en:<12}{step_ko[step]:<9}{zn_kegg:>10}{zn_og:>9}{nip_og:>8}")
    rows.append((en,step_ko[step],zn_kegg,zn_og,nip_og)); t1+=zn_kegg; t2+=zn_og; t3+=nip_og
print(f"{'TOTAL':<12}{'':<9}{t1:>10}{t2:>9}{t3:>8}")
with open("pathway_copynumber_ZN65_vs_NIP.tsv","w") as fh:
    fh.write("Enzyme_step\tKO\tZN65_KEGG_copies\tZN65_orthogroup_copies\tNipponbare_orthogroup_copies\n")
    for r in rows: fh.write("\t".join(str(x) for x in r)+"\n")
    fh.write(f"TOTAL\t\t{t1}\t{t2}\t{t3}\n")
# correlation-ish summary
import statistics
ratios=[ (r[3]/r[4]) for r in rows if r[4]>0 ]
print(f"\nZN65_OG/NIP_OG ratio: median={statistics.median(ratios):.2f}, mean={statistics.mean(ratios):.2f} (1.0 = conserved)")
print("Wrote pathway_copynumber_ZN65_vs_NIP.tsv")
