source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
A12g="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/A12_genetree2"
GEN="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/genomes"
B="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定"
OLD="/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/A12_genetree"
mkdir -p "$A12g/genes"; cd "$A12g"
CDS="/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.cds"
echo "== gene-tree2 START $(date) =="
# more candidates: every 30th gene, length 1000-2500
seqkit fx2tab -nl "$CDS" 2>/dev/null | awk '$2>=1000 && $2<=2500{print $1}' | awk 'NR%30==0' > cand.ids
seqkit grep -f cand.ids "$CDS" 2>/dev/null > cand.cds
echo "candidates: $(grep -c '>' cand.cds)"
# reuse blastdbs from genetree run
for nm in ZN65 Nipponbare MH63 ZS97 CempoIreng N22 Orufipogon; do
  [ -f "$OLD/db_$nm.nsq" ] && ln -sf "$OLD/db_$nm".* . 2>/dev/null
done
# batch blastn: one call per genome
for nm in ZN65 Nipponbare MH63 ZS97 CempoIreng N22 Orufipogon; do
  blastn -query cand.cds -db db_$nm -outfmt '6 qseqid sseqid sstart send pident length qlen bitscore' -evalue 1e-30 -max_target_seqs 3 -num_threads 8 > hits_$nm.tsv 2>/dev/null
  echo "blast $nm done: $(wc -l < hits_$nm.tsv) hits"
done
python3 - <<'PY'
import os,subprocess
panel=["ZN65","Nipponbare","MH63","ZS97","CempoIreng","N22","Orufipogon"]
gpath={"ZN65":"/mnt/data2/墨江紫米研究/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa",
"Nipponbare":"/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/08_全基因组SV_ZN65vsNIP/nip_genome.fa",
"MH63":"/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/10_多基因组Kala4/MH63RS3.fasta",
"ZS97":"/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/10_多基因组Kala4/ZS97RS3.fasta",
"CempoIreng":"/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/genomes/CempoIreng_black.fna",
"N22":"/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/genomes/N22_aus.fna",
"Orufipogon":"/mnt/data2/墨江紫米研究/07_分析_Os02g基因鉴定/22_泛基因组_A10A12/genomes/Orufipogon_wild.fna"}
def besthit(nm):
    d={}
    for ln in open(f"hits_{nm}.tsv"):
        q,s,ss,se,pid,ln_,qlen,bit=ln.split("\t")
        pid=float(pid); aln=int(ln_); qlen=int(qlen); bit=float(bit)
        if pid<88 or aln<0.7*qlen: continue
        d.setdefault(q,[]).append((bit,s,int(ss),int(se)))
    res={}
    for q,h in d.items():
        h.sort(reverse=True)
        if len(h)>=2 and h[1][0]>0.85*h[0][0]: continue  # not single-copy
        res[q]=h[0]
    return res
hits={nm:besthit(nm) for nm in panel}
common=set(hits["ZN65"])
for nm in panel: common&=set(hits[nm])
print("single-copy genes in ALL 7:",len(common))
comp=str.maketrans("ACGTacgtNn","TGCAtgcaNn")
kept=0
for g in sorted(common):
    recs=[]
    ok=True
    for nm in panel:
        bit,s,ss,se=hits[nm][g]
        a,b=min(ss,se),max(ss,se); rev=ss>se
        sub=subprocess.run(["samtools","faidx",gpath[nm],f"{s}:{a}-{b}"],capture_output=True,text=True).stdout
        sq="".join(sub.split("\n")[1:])
        if not sq: ok=False; break
        if rev: sq=sq.translate(comp)[::-1]
        recs.append((nm,sq))
    if not ok or len(recs)!=7: continue
    with open(f"genes/{g}.fa","w") as o:
        for nm,sq in recs: o.write(f">{nm}\n{sq}\n")
    kept+=1
print("genes extracted:",kept)
PY
# align + concat
> concat.list
for f in genes/*.fa; do
  [ $(grep -c '>' "$f") -eq 7 ] || continue
  mafft --auto --thread 8 "$f" > "${f%.fa}.aln" 2>/dev/null && echo "${f%.fa}.aln" >> concat.list
done
echo "genes aligned: $(wc -l < concat.list)"
python3 - <<'PY'
order=["ZN65","Nipponbare","MH63","ZS97","CempoIreng","N22","Orufipogon"]
cat={s:[] for s in order}
for aln in open("concat.list"):
    aln=aln.strip(); seqs={}; cur=None
    for ln in open(aln):
        if ln[0]==">": cur=ln[1:].strip(); seqs[cur]=[]
        else: seqs[cur].append(ln.strip())
    seqs={k:"".join(v) for k,v in seqs.items()}
    L=len(next(iter(seqs.values())))
    for s in order: cat[s].append(seqs.get(s,"-"*L))
open("concat.fasta","w").write("".join(f">{s}\n{''.join(cat[s])}\n" for s in order))
print("concat length:",len("".join(cat['ZN65'])))
PY
iqtree -s concat.fasta -m GTR+G -bb 1000 -nt 8 -redo -pre gene_tree2 -o Orufipogon 2>/dev/null
echo "=== GENE TREE2 (outgroup=Orufipogon) ==="; cat gene_tree2.contree
echo "== DONE $(date) =="; touch "$A12g/DONE2"
