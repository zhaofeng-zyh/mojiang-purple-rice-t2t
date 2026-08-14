# ⛔ SUPERSEDED —— 本脚本会按**已被推翻的口径**再生产物，入口硬失败（2026-07-29 加装）
#
# 已被推翻的旧口径：「OsB2 启动子有一段 ~5.9 kb / ~6.3 kb 的 ZN65 特异插入，白等位缺失，
# 可用存在-缺失共显性分型」。现行口径：OsB2 **上游**是 ~5.87 kb block A 的 LINE 关联
# **复合局部片段重复** —— ZN65 3 拷贝（A1/A2/A3），日本晴/MH63/ZS97/N22/Cempo Ireng/
# 野生 O. rufipogon **6 个对照各单拷贝**（99.3–100%）；该片段**在所有 7 个受检基因组中都存在**，
# 差异只是拷贝数 3 vs 1，**不存在「对照缺失」**，故**不能做存在-缺失共显性分型**。
# 接头缺 TSD 与完整末端重复 → **不是逆转座子插入**。内含子 LINE1-11_OS 实测 **5,934 bp**（非 6.3 kb）；
# 三拷贝区 A1 起–A3 止 **115,450 bp**（非 ~119 kb）。
#
# 保留本文件仅为留痕。若确需运行（例如复现旧图以作对照），显式设
# `ALLOW_SUPERSEDED=1` —— 但产物**不得**作为现行结论使用。
import os as _os, sys as _sys
if _os.environ.get("ALLOW_SUPERSEDED") != "1":
    _sys.exit("⛔ 本脚本按已被推翻的口径产出（详见文件头注释）。"
              "现行口径见 12_论文Paper1_Manuscript/_投稿包_TheCropJournal_20260724/。"
              "确需留痕重现请设 ALLOW_SUPERSEDED=1。")
import sys
M="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/manuscript.tex"
s=open(M,encoding="utf-8").read()
old=("Orthologues were identified by reciprocal BLAST against\n"
"the IRGSP-1.0 proteome.")
new=("Orthologues were identified by reciprocal BLAST against\n"
"the IRGSP-1.0 proteome. Presence or absence of the ZN65-specific \\textit{OsB2}\n"
"distal-promoter insertion was assessed across a seven-genome panel that additionally\n"
"included N22 (aus; GenBank GCA\\_001952365.2) and \\textit{Oryza rufipogon} W1943\n"
"(wild progenitor; GCA\\_000817225.1); the orthologous \\textit{OsB2} locus was located\n"
"in each genome by anchoring on the 34-kb conserved region, and the 5.9-kb insertion\n"
"was scored only within the orthologous interval to avoid false positives from\n"
"dispersed transposon copies.")
n=s.count(old)
if n!=1: print("ERROR count",n); sys.exit(1)
open(M,"w",encoding="utf-8").write(s.replace(old,new))
print("Methods: A10 panel + orthology-anchored test added OK")
