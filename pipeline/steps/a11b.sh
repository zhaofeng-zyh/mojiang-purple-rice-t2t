#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate cgsv
P=/mnt/data2/墨江紫米研究
PEP=$P/01_基因组_组装与注释_Genome/02.annotation/02.gene/ZN65.longest.pep
LIST=$P/07_分析_Os02g基因鉴定/09_花青素通路清单/ZN65_花青素通路基因清单.tsv
W=$P/07_分析_Os02g基因鉴定/18_通路拷贝数_A11; cd $W
seqkit grep -r -f genes.txt $PEP > zn65_pathway.pep 2>/dev/null
echo "proteins extracted: $(grep -c '>' zn65_pathway.pep)"
blastp -query zn65_pathway.pep -db nipdb -evalue 1e-5 -max_target_seqs 1 -outfmt "6 qseqid sseqid pident length qlen" -num_threads 8 2>/dev/null | sort -k1,1 -k3,3nr | awk "!seen[\$1]++" > best.tsv
echo "blast hits: $(wc -l < best.tsv)"
python3 - <<PY
step={}
for ln in open("$LIST"):
    p=ln.rstrip("\n").split("\t")
    if p[0]=="通路步骤": continue
    step[p[2]]=p[0]
best={}
for ln in open("best.tsv"):
    f=ln.rstrip("\n").split("\t"); g=f[0].split(".")[0]; best[g]=float(f[2])
from collections import defaultdict
zc=defaultdict(int); orth=defaultdict(int); ids=defaultdict(list)
for g,st in step.items():
    zc[st]+=1
    b=best.get(g) or best.get(g+".1")
    # match by prefix
    if b is None:
        for k,v in best.items():
            if k.startswith(g): b=v; break
    if b and b>=40: orth[st]+=1; ids[st].append(b)
print("%-30s %-5s %-8s %s"%("Pathway step","ZN65","NIP_orth","mean_id%"))
for st in sorted(zc, key=lambda x: -zc[x]):
    mi=sum(ids[st])/len(ids[st]) if ids[st] else 0
    en=st.split()[0]
    print("%-30s %-5d %-8d %.1f"%(en, zc[st], orth[st], mi))
print("TOTAL:", sum(zc.values()), " w/ NIP ortholog:", sum(orth.values()))
PY
