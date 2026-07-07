source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
A12g="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/A12_genetree"
GEN="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/genomes"
B="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定"
mkdir -p "$A12g/genes"; cd "$A12g"
ZN65="/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa"
NIP="$B/08_全基因组SV_ZN65vsNIP/nip_genome.fa"
MH63="$B/10_多基因组Kala4/MH63RS3.fasta"; ZS97="$B/10_多基因组Kala4/ZS97RS3.fasta"
CDS="/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.cds"
echo "== gene-tree START $(date) =="
# 1) candidate conserved single-copy CDS: deterministic sample, length 1200-2400
seqkit fx2tab -nl "$CDS" 2>/dev/null | awk '$2>=1200 && $2<=2400{print $1}' | awk 'NR%150==0' | head -90 > cand.ids
seqkit grep -f cand.ids "$CDS" 2>/dev/null > cand.cds
echo "candidate genes: $(grep -c '>' cand.cds)"
# 2) blastdb for each genome
declare_panel="ZN65:$ZN65 Nipponbare:$NIP MH63:$MH63 ZS97:$ZS97 CempoIreng:$GEN/CempoIreng_black.fna N22:$GEN/N22_aus.fna Orufipogon:$GEN/Orufipogon_wild.fna"
for e in $declare_panel; do nm=${e%%:*}; g=${e#*:}; [ -f db_$nm.nsq ] || makeblastdb -in "$g" -dbtype nucl -out db_$nm >/dev/null 2>&1; done
echo "dbs ready $(date)"
# 3) for each candidate gene, get single-copy strong hit in each genome -> extract
python3 - <<'PY'
import subprocess,os
panel=[("ZN65","db_ZN65"),("Nipponbare","db_Nipponbare"),("MH63","db_MH63"),("ZS97","db_ZS97"),("CempoIreng","db_CempoIreng"),("N22","db_N22"),("Orufipogon","db_Orufipogon")]
gpath={"ZN65":"/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa",
"Nipponbare":"/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP/nip_genome.fa",
"MH63":"/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/10_多基因组Kala4/MH63RS3.fasta",
"ZS97":"/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/10_多基因组Kala4/ZS97RS3.fasta",
"CempoIreng":"/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/genomes/CempoIreng_black.fna",
"N22":"/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/genomes/N22_aus.fna",
"Orufipogon":"/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/genomes/Orufipogon_wild.fna"}
# read candidate cds
genes={}; cur=None
for ln in open("cand.cds"):
    if ln[0]==">": cur=ln[1:].split()[0]; genes[cur]=[]
    else: genes[cur].append(ln.strip())
genes={k:"".join(v) for k,v in genes.items()}
open("g.fa","w")
kept=0
for gid,seq in genes.items():
    open("g.fa","w").write(f">{gid}\n{seq}\n")
    ok=True; recs={}
    for nm,db in panel:
        out=subprocess.run(["blastn","-query","g.fa","-db",db,"-outfmt","6 sseqid sstart send pident length qlen bitscore","-evalue","1e-30","-max_target_seqs","5"],capture_output=True,text=True).stdout.strip().split("\n")
        hits=[l.split("\t") for l in out if l]
        if not hits: ok=False; break
        hits=[h for h in hits if float(h[3])>=85 and int(h[4])>=0.8*int(h[5])]
        if not hits: ok=False; break
        hits.sort(key=lambda h:-float(h[6]))
        if len(hits)>=2 and float(hits[1][6])>0.8*float(hits[0][6]): ok=False; break  # not single-copy
        sid,ss,se=hits[0][0],int(hits[0][1]),int(hits[0][2])
        recs[nm]=(sid,min(ss,se),max(ss,se),ss>se)
    if not ok: continue
    # extract & write per-gene multifasta
    with open(f"genes/{gid}.fa","w") as o:
        for nm,db in panel:
            sid,a,b,rev=recs[nm]
            sub=subprocess.run(["samtools","faidx",gpath[nm],f"{sid}:{a}-{b}"],capture_output=True,text=True).stdout
            sq="".join(sub.split("\n")[1:])
            if rev:
                comp={"A":"T","T":"A","G":"C","C":"G","N":"N","a":"t","t":"a","g":"c","c":"g","n":"n"}
                sq="".join(comp.get(c,"N") for c in reversed(sq))
            o.write(f">{nm}\n{sq}\n")
    kept+=1
print("single-copy genes kept:",kept)
PY
# 4) align each gene, keep those with all 7 taxa, concatenate
> concat.list
for f in genes/*.fa; do
  [ $(grep -c '>' "$f") -eq 7 ] || continue
  mafft --auto --thread 8 "$f" > "${f%.fa}.aln" 2>/dev/null && echo "${f%.fa}.aln" >> concat.list
done
echo "genes aligned (7 taxa): $(wc -l < concat.list)"
# concatenate with seqkit
python3 - <<'PY'
import glob
order=["ZN65","Nipponbare","MH63","ZS97","CempoIreng","N22","Orufipogon"]
cat={s:[] for s in order}
for aln in open("concat.list"):
    aln=aln.strip()
    seqs={}; cur=None
    for ln in open(aln):
        if ln[0]==">": cur=ln[1:].strip(); seqs[cur]=[]
        else: seqs[cur].append(ln.strip())
    seqs={k:"".join(v) for k,v in seqs.items()}
    L=len(next(iter(seqs.values())))
    for s in order: cat[s].append(seqs.get(s,"-"*L))
open("concat.fasta","w").write("".join(f">{s}\n{''.join(cat[s])}\n" for s in order))
print("concat length:", len("".join(cat['ZN65'])))
PY
iqtree -s concat.fasta -m GTR+G -bb 1000 -nt 8 -redo -pre gene_tree 2>/dev/null
echo "=== GENE TREE ==="; cat gene_tree.contree
echo "== gene-tree DONE $(date) =="
touch "$A12g/GENETREE_DONE"
