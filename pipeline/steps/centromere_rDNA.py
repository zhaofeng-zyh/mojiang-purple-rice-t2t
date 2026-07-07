#!/usr/bin/env python3
"""Localize centromeres (CentO satellite arrays) and 45S/5S rDNA arrays in the ZN65 T2T assembly.
Inputs (existing project files); outputs tables into the submission analysis dir."""
import re, collections, os
BASE="/mnt/data2/墨江紫米研究"
LAND=f"{BASE}/07_分析_Os02g基因鉴定/11_基因组景观"
TRF=f"{BASE}/01_基因组_组装与注释_Genome/02.annotation/01.repeat/trf.gff"
RRNA=f"{BASE}/01_基因组_组装与注释_Genome/02.annotation/03.ncRNA/rRNA.gff"
OUT=f"{BASE}/12_论文Paper1_Manuscript/_SUBMISSION_PlantComm/analysis"
CHRS=[f"Chr{i}" for i in range(1,13)]
sizes={}
for l in open(f"{LAND}/genome.txt"):
    c,s=l.split(); sizes[c]=int(s)

# ---------- 1. Centromere via CentO satellite (rice CentO ~155-165 bp) ----------
# collect tandem-repeat spans with period 150-170 bp; cluster per chr; centromere = densest 1-Mb region
cento=collections.defaultdict(list)   # chr -> list of (start,end,span_bp,period)
allsat=collections.defaultdict(list)
for l in open(TRF):
    if l.startswith('#'): continue
    f=l.rstrip('\n').split('\t')
    if len(f)<9: continue
    c=f[0]; 
    if c not in sizes: continue
    st,en=int(f[3]),int(f[4])
    m=re.search(r'PeriodSize=(\d+)',f[8]); 
    if not m: continue
    per=int(m.group(1))
    span=en-st
    allsat[c].append((st,en,span,per))
    if 150<=per<=170 and span>=1000:
        cento[c].append((st,en,span,per))

def densest_region(spans, win=1_000_000):
    """return (start,end,total_satellite_bp) of the win-sized region with max satellite bp"""
    if not spans: return None
    # bin by 100kb, sum span, then sliding 1Mb
    binsz=100_000
    b=collections.Counter()
    for st,en,sp,per in spans:
        b[st//binsz]+=sp
    if not b: return None
    maxk=max(b); 
    best=(0,0,-1)
    for start in range(0, (maxk+1)*binsz, binsz):
        tot=sum(b.get(k,0) for k in range(start//binsz,(start+win)//binsz))
        if tot>best[2]: best=(start,start+win,tot)
    return best

print("=== CENTROMERE (CentO 155-165bp satellite arrays) ===")
cen_table=[]
for c in CHRS:
    reg=densest_region(cento[c])
    if reg is None or reg[2]<5000:
        # fallback: max TE-coverage 500kb window
        reg2=None
    # also compute CentO monomer consensus most common period
    pers=[p for (_,_,_,p) in cento[c]]
    core = densest_region(cento[c])
    tot_cento=sum(sp for (_,_,sp,_) in cento[c])
    if core:
        st,en,bp=core
        # clip to chr
        en=min(en,sizes[c])
        cen_table.append((c,st,en,(st+en)//2,tot_cento,bp))
        print(f"{c}\tcentromere~{st:,}-{en:,}\tmid={(st+en)//2:,}\tCentO_bp_in_core={bp:,}\tCentO_total={tot_cento:,}")
    else:
        cen_table.append((c,'NA','NA','NA',tot_cento,0))
        print(f"{c}\tno strong CentO array (total CentO bp={tot_cento})")

# ---------- 2. rDNA arrays (45S: 18S/5.8S/25S tandem; 5S separate) ----------
rrna=collections.defaultdict(list)   # chr -> (start,end,type)
typ_count=collections.Counter()
for l in open(RRNA):
    if l.startswith('#'): continue
    f=l.rstrip('\n').split('\t')
    if len(f)<9 or f[0] not in sizes: continue
    st,en=int(f[3]),int(f[4])
    ann=f[8]
    # classify by rRNA subunit
    if '5.8S' in ann or '5_8S' in ann: t='45S(5.8S)'
    elif '18S' in ann or 'SSU' in ann: t='45S(18S)'
    elif '25S' in ann or '28S' in ann or '26S' in ann or 'LSU' in ann: t='45S(25S)'
    elif '5S' in ann: t='5S'
    else: t='other'
    rrna[f[0]].append((st,en,t)); typ_count[t]+=1
print("\n=== rRNA feature types ===")
for t,n in typ_count.most_common(): print(f"  {t}: {n}")

def cluster(items, gap=50000):
    items=sorted(items)
    clusters=[]; cur=[items[0]]
    for it in items[1:]:
        if it[0]-cur[-1][1] <= gap: cur.append(it)
        else: clusters.append(cur); cur=[it]
    clusters.append(cur)
    return clusters

print("\n=== rDNA ARRAYS (features clustered, gap<=50kb) ===")
rdna_table=[]
for c in CHRS:
    if c not in rrna: continue
    # split 45S vs 5S
    for label,pred in [('45S', lambda t:t.startswith('45S')),('5S', lambda t:t=='5S')]:
        pts=[(s,e) for (s,e,t) in rrna[c] if pred(t)]
        if not pts: continue
        for cl in cluster(pts):
            st=min(p[0] for p in cl); en=max(p[1] for p in cl); n=len(cl)
            if (en-st)>=20000 and n>=10:  # substantial array
                rdna_table.append((c,label,st,en,en-st,n))
                print(f"{c}\t{label}\t{st:,}-{en:,}\t{(en-st)/1000:.0f}kb\tunits={n}")

# ---------- write tables ----------
with open(f"{OUT}/centromere_table.tsv","w") as fh:
    fh.write("Chromosome\tCentromere_start\tCentromere_end\tCentromere_mid\tCentO_bp_total\tCentO_bp_core1Mb\n")
    for r in cen_table: fh.write("\t".join(str(x) for x in r)+"\n")
with open(f"{OUT}/rDNA_table.tsv","w") as fh:
    fh.write("Chromosome\trDNA_class\tArray_start\tArray_end\tArray_len_bp\tN_units\n")
    for r in rdna_table: fh.write("\t".join(str(x) for x in r)+"\n")
print(f"\nWrote centromere_table.tsv ({len(cen_table)} chr) and rDNA_table.tsv ({len(rdna_table)} arrays)")
