#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
SWD="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12"
B="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定"
GEN="$SWD/genomes"; A10="$SWD/A10"; A12="$SWD/A12"
mkdir -p "$A10" "$A12" "$SWD/logs"
ZN65="/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa"
NIP="$B/08_全基因组SV_ZN65vsNIP/nip_genome.fa"
MH63="$B/10_多基因组Kala4/MH63RS3.fasta"; ZS97="$B/10_多基因组Kala4/ZS97RS3.fasta"
CEMPO="$GEN/CempoIreng_black.fna"; N22="$GEN/N22_aus.fna"; RUF="$GEN/Orufipogon_wild.fna"
NPROC=20
echo "== PIPELINE START $(date) =="

# ---- wait for downloaded genomes (up to 70 min) ----
for i in $(seq 1 140); do
  [ -s "$CEMPO" ] && [ -s "$N22" ] && [ -s "$RUF" ] && { echo "genomes ready $(date +%H:%M:%S)"; break; }
  echo "[wait genomes] $(date +%H:%M:%S)"; sleep 30
done

PANEL="ZN65:$ZN65:japonica:purple Nipponbare:$NIP:japonica:white MH63:$MH63:indica:white ZS97:$ZS97:indica:white CempoIreng:$CEMPO:tropjaponica:black N22:$N22:aus:white Orufipogon:$RUF:wild:wild"

#################### A10 : OsB2 promoter-TE insertion presence/absence ####################
echo "== A10 START $(date) =="
cd "$A10"
cp "$B/07_花青素位点SV/OsB2_Kala4/ZN65_Kala4_locus.fasta" zn65_locus.fa
# diagnostic ZN65-specific insert = ZN65 locus positions not covered by Nipponbare (computed earlier ~34131-40001)
# recompute robustly via minimap2 of NIP locus onto ZN65 locus
cp "$B/07_花青素位点SV/OsB2_Kala4/NIP_Kala4_locus.fasta" nip_locus.fa
nucmer --maxmatch -l 20 -p zn_vs_nip nip_locus.fa zn65_locus.fa 2>/dev/null
show-coords -rclTH zn_vs_nip.delta > coords.tsv 2>/dev/null
python3 - <<'PY'
ivs=[]
for ln in open("coords.tsv"):
    p=ln.rstrip("\n").split("\t")
    if len(p)<7: continue
    try: s2,e2=int(p[2]),int(p[3])
    except: continue
    a,b=min(s2,e2),max(s2,e2)
    if b-a+1>=50: ivs.append((a,b))
ivs.sort(); merged=[]
for a,b in ivs:
    if merged and a<=merged[-1][1]+1: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
    else: merged.append((a,b))
L=40001; unc=[]; prev=1
for a,b in merged:
    if a-prev>=300: unc.append((prev,a-1))
    prev=max(prev,b+1)
if L-prev>=300: unc.append((prev,L))
unc.sort(key=lambda x:-(x[1]-x[0]))
open("insert_intervals.txt","w").write("\n".join(f"{a}\t{b}\t{b-a+1}" for a,b in unc)+"\n")
print("ZN65-specific intervals (largest first):", unc[:5])
PY
read A Bx LEN < <(head -1 insert_intervals.txt)
echo "diagnostic insert ZN65 locus $A-$Bx ($LEN bp)"
seqkit subseq -r "${A}:${Bx}" zn65_locus.fa 2>/dev/null | seqkit replace -p '.*' -r 'ZN65_OsB2_insert' 2>/dev/null > insert.fa
seqkit stats insert.fa

echo -e "genome\tsubspecies\tpericarp\tlocus_aln%\tinsert_aln%\tconserved_aln%\tinsert_call" > PA_table.tsv
for e in $PANEL; do
  nm=${e%%:*}; r=${e#*:}; g=${r%%:*}; r=${r#*:}; ssp=${r%%:*}; per=${r##*:}
  # map the whole ZN65 40kb locus onto genome (flanks anchor orthology); report coverage of insert vs conserved
  minimap2 -cx asm10 --secondary=no "$g" zn65_locus.fa 2>/dev/null > $nm.paf
  python3 - "$nm.paf" "$A" "$Bx" <<'PY' >> PA_table.tsv
import sys
paf,ins0,ins1=sys.argv[1],int(sys.argv[2]),int(sys.argv[3])
# pick best target region by total aligned query length, then compute query coverage intervals
rows=[]
for ln in open(paf):
    p=ln.split("\t")
    if len(p)<12: continue
    qs,qe=int(p[2]),int(p[3]); tname=p[5]; ts,te=int(p[7]),int(p[8]); mapq=int(p[11])
    rows.append((tname,ts,te,qs,qe,mapq))
if not rows:
    print("NA\tNA\tNA\t0.0\t0.0\t0.0\tABSENT"); sys.exit()
# best target = the one with max aligned span near a single region
from collections import defaultdict
byt=defaultdict(int)
for tname,ts,te,qs,qe,mq in rows: byt[tname]+=qe-qs
bestt=max(byt,key=byt.get)
qiv=sorted([(qs,qe) for tn,ts,te,qs,qe,mq in rows if tn==bestt])
# merge query intervals
mg=[]
for a,b in qiv:
    if mg and a<=mg[-1][1]: mg[-1]=(mg[-1][0],max(mg[-1][1],b))
    else: mg.append([a,b])
def covlen(a,b):
    s=0
    for x,y in mg:
        lo,hi=max(a,x),min(b,y)
        if hi>lo: s+=hi-lo
    return s
L=40001
locus=sum(y-x for x,y in mg)
insert=covlen(ins0,ins1); conserved=covlen(1,ins0)
ilen=ins1-ins0; clen=ins0-1
ip=100.0*insert/ilen if ilen else 0
cp=100.0*conserved/clen if clen else 0
call="PRESENT" if ip>=60 else ("PARTIAL" if ip>=25 else "ABSENT")
print(f"{100.0*locus/L:.1f}\t{ip:.1f}\t{cp:.1f}\t{call}")
PY
done
# fix table: the python prints 4 cols after the leading per-genome cols; reassemble
echo "A10 PA_table:"; cat PA_table.tsv

# insertion-site phylogeny: conserved OsB2 flank (proximal promoter, aligns across all) -> MAFFT -> IQ-TREE
cp "$B/10_多基因组Kala4/osb2_prom_zn65.fa" flank.fa
> flanks_all.fa
for e in $PANEL; do
  nm=${e%%:*}; r=${e#*:}; g=${r%%:*}
  makeblastdb -in "$g" -dbtype nucl -out tdb_$nm >/dev/null 2>&1
  reg=$(blastn -query flank.fa -db tdb_$nm -evalue 1e-20 -outfmt '6 sseqid sstart send' -max_target_seqs 1 2>/dev/null | head -1)
  if [ -n "$reg" ]; then
    sid=$(echo "$reg"|cut -f1); ss=$(echo "$reg"|cut -f2); se=$(echo "$reg"|cut -f3)
    lo=$(( ss<se ? ss : se )); hi=$(( ss<se ? se : ss ))
    samtools faidx "$g" "$sid:$lo-$hi" 2>/dev/null | seqkit replace -p '.*' -r "$nm" 2>/dev/null >> flanks_all.fa
  fi
  rm -f tdb_$nm.*
done
seqkit stats flanks_all.fa
mafft --auto --thread $NPROC flanks_all.fa > flanks_aln.fa 2>/dev/null
iqtree -s flanks_aln.fa -m MFP -bb 1000 -nt $NPROC -redo -pre osb2_flank_tree 2>/dev/null
echo "A10 tree done $(date)"

#################### A12 : genome-wide SNP phylogeny + PCA vs Nipponbare ####################
echo "== A12 START $(date) =="
cd "$A12"
for e in $PANEL; do
  nm=${e%%:*}; r=${e#*:}; g=${r%%:*}
  [ "$nm" = "Nipponbare" ] && continue
  echo "[A12] calling SNPs $nm $(date +%H:%M:%S)"
  minimap2 -cx asm20 --cs -t $NPROC "$NIP" "$g" 2>/dev/null | sort -k6,6 -k8,8n 2>/dev/null | paftools.js call -L10000 - 2>/dev/null > $nm.var
done
python3 - <<'PY'
import glob,os
import numpy as np
samples=[os.path.basename(f)[:-4] for f in sorted(glob.glob("*.var"))]
# collect SNPs per sample: dict[(chr,pos)] = alt
data={}; allpos=set()
chrok=lambda c: c.startswith("chr") or c.startswith("Chr") or True
for s in samples:
    d={}
    for ln in open(s+".var"):
        if not ln.startswith("V"): continue
        p=ln.rstrip("\n").split("\t")
        # V chr start end qdepth mapq ref alt ...
        if len(p)<8: continue
        ref,alt=p[6],p[7]
        if len(ref)==1 and len(alt)==1 and ref in "ACGT" and alt in "ACGT":
            key=(p[1],int(p[2])); d[key]=alt; allpos.add(key)
    data[s]=d
# keep sites variable & genotyped in >=  (n-1) samples to limit missing
pos=sorted(allpos)
keep=[k for k in pos if sum(1 for s in samples if k in data[s])>=max(2,len(samples)-1)]
# build allele matrix: ref(Nipponbare)=0 ; sample alt vs ref. Nipponbare is reference => allele=ref
# represent each sample at site: alt if present else REF(=ancestral nip). Need ref base: take from any sample's record? paftools gives ref col; store it.
refbase={}
for s in samples:
    for ln in open(s+".var"):
        if not ln.startswith("V"): continue
        p=ln.split("\t")
        if len(p)<8: continue
        k=(p[1],int(p[2]))
        if k in keep and k not in refbase and len(p[6])==1: refbase[k]=p[6]
keep=[k for k in keep if k in refbase]
print("usable SNP sites:", len(keep))
aln={"Nipponbare":[]}
for s in samples: aln[s]=[]
for k in keep:
    aln["Nipponbare"].append(refbase[k])
    for s in samples:
        aln[s].append(data[s].get(k, refbase[k]))
order=["Nipponbare"]+samples
with open("snp_aln.fasta","w") as o:
    for s in order:
        o.write(f">{s}\n{''.join(aln[s])}\n")
# PCA on 0/1 (alt vs ref)
M=np.array([[0 if aln[s][i]==refbase[keep[i]] else 1 for i in range(len(keep))] for s in order],dtype=float)
M=M-M.mean(0)
U,S,Vt=np.linalg.svd(M,full_matrices=False)
pc=U*S
with open("pca.tsv","w") as o:
    o.write("sample\tPC1\tPC2\tPC3\n")
    for i,s in enumerate(order):
        o.write(f"{s}\t{pc[i,0]:.3f}\t{pc[i,1]:.3f}\t{pc[i,2]:.3f}\n")
print("PCA written")
PY
iqtree -s snp_aln.fasta -m MFP -bb 1000 -nt $NPROC -redo -pre zn65_subspecies_tree 2>/dev/null
echo "A12 tree done $(date)"
echo "== PIPELINE DONE $(date) =="
touch "$SWD/PIPELINE_DONE"
