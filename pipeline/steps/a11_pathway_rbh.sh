#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
echo "A11 START: $(date '+%Y-%m-%d %H:%M:%S')"
P=/mnt/data2/墨江紫米研究
PEP=$P/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.pep
NIP=$P/07_分析_Os02g基因鉴定/02_参考序列/nip_pep.fa
LIST=$P/07_分析_Os02g基因鉴定/09_花青素通路清单/ZN65_花青素通路基因清单.tsv
W=$P/07_分析_Os02g基因鉴定/18_通路拷贝数_A11; mkdir -p $W; cd $W
# extract ZN65 pathway proteins
tail -n+2 "$LIST" | cut -f3 > genes.txt
seqkit grep -f genes.txt $PEP > zn65_pathway.pep 2>/dev/null
echo "pathway proteins extracted: $(grep -c '>' zn65_pathway.pep)"
# blast db of Nipponbare
makeblastdb -in $NIP -dbtype prot -out nipdb > /dev/null 2>&1
# blastp pathway -> Nipponbare, best hit
blastp -query zn65_pathway.pep -db nipdb -evalue 1e-5 -max_target_seqs 1 -outfmt "6 qseqid sseqid pident length qlen" -num_threads 8 2>/dev/null | sort -k1,1 -k3,3nr | awk "!seen[\$1]++" > best.tsv
# join enzyme step + count
python3 - <<PY
step={}; 
for ln in open("$LIST"):
    p=ln.rstrip("\n").split("\t")
    if p[0]=="通路步骤": continue
    step[p[2]]=p[0]
best={}
for ln in open("best.tsv"):
    q,s,pid,length,qlen=ln.rstrip("\n").split("\t")
    best[q]=(s,float(pid))
from collections import defaultdict
zc=defaultdict(int); orth=defaultdict(int); ids=defaultdict(list)
for g,st in step.items():
    zc[st]+=1
    if g in best and best[g][1]>=40:
        orth[st]+=1; ids[st].append(best[g][1])
print("%-34s %-6s %-10s %s"%("Pathway step","ZN65","NIP_orth","mean_id%"))
for st in sorted(zc, key=lambda x: -zc[x]):
    mi=sum(ids[st])/len(ids[st]) if ids[st] else 0
    print("%-34s %-6d %-10d %.1f"%(st[:33], zc[st], orth[st], mi))
print("TOTAL ZN65 pathway genes:", sum(zc.values()), " with NIP ortholog(>=40% id):", sum(orth.values()))
PY
echo "A11 FINISH: $(date '+%Y-%m-%d %H:%M:%S')"
echo "A11_DONE"
