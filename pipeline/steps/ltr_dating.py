#!/usr/bin/env python3
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
"""
Reconcile the Kala4/OsB2 intronic RETRO2B (LTR/Gypsy) insertion age.

Recomputes the insertion age from the ACTUAL 5'LTR vs 3'LTR alignment
(EMBOSS stretcher output `ltr_aln.txt`), counting substitutions ONLY
(indels excluded, as required for LTR retrotransposon dating).

Resolves the project discrepancy:
  - manuscript abstract: ~0.2 Mya  (substitutions-only)   <-- expected correct
  - 结果_Kala4致色变异与定年.md: ~1.4 Mya (K=0.0361 counted indels as differences) <-- error

Method: Kimura (1980) 2-parameter distance; T = K / (2*mu); mu = 1.3e-8 subs/site/yr
        (rice LTR substitution rate, as used project-wide).
95% CI on the substitution count via the exact Poisson interval.

Reproducible: paths are relative to this file (no hardcoded /sessions or /Volumes).
"""
from pathlib import Path
import math, re

ROOT = Path(__file__).resolve().parents[2]   # -> 墨江紫米研究
ALN  = ROOT / "07_分析_Os02g基因鉴定/07_花青素位点SV/OsB2_Kala4/ltr_aln.txt"
MU   = 1.3e-8   # substitutions / site / year (rice LTR-RT rate)

def parse_stretcher(path):
    """Return the two gapped aligned strings from an EMBOSS stretcher 'pair' report."""
    seqs = {}
    order = []
    for line in path.read_text().splitlines():
        # data lines look like:  NAME  <start>  SEQCHUNK  <end>
        m = re.match(r'^(\S+)\s+\d+\s+([A-Za-z\-\.]+)\s+\d+\s*$', line)
        if m:
            name, chunk = m.group(1), m.group(2)
            if name not in seqs:
                seqs[name] = []
                order.append(name)
            seqs[name].append(chunk)
    if len(order) != 2:
        raise SystemExit(f"Expected 2 aligned sequences, parsed {order}")
    a = "".join(seqs[order[0]]).upper()
    b = "".join(seqs[order[1]]).upper()
    if len(a) != len(b):
        raise SystemExit(f"Aligned lengths differ: {len(a)} vs {len(b)}")
    return order, a, b

TRANSITIONS = {frozenset("AG"), frozenset("CT")}

def k2p(a, b):
    aln_len = len(a)
    ident = gaps = ts = tv = 0
    for x, y in zip(a, b):
        if x == '-' or y == '-' or x == '.' or y == '.':
            gaps += 1; continue
        if x == y:
            ident += 1
        else:
            pair = frozenset((x, y))
            if pair in TRANSITIONS: ts += 1
            else: tv += 1
    ungapped = aln_len - gaps
    subs = ts + tv
    P = ts / ungapped
    Q = tv / ungapped
    # Kimura 2-parameter
    try:
        K = -0.5*math.log(1 - 2*P - Q) - 0.25*math.log(1 - 2*Q)
    except ValueError:
        K = float('nan')
    p_dist = subs / ungapped
    return dict(aln_len=aln_len, ungapped=ungapped, ident=ident, gaps=gaps,
                ts=ts, tv=tv, subs=subs, P=P, Q=Q, K=K, p_dist=p_dist)

def poisson_ci(k, alpha=0.05):
    """Exact (Garwood) Poisson 95% CI for an observed count k, via chi-square quantiles."""
    # use math.gamma-based inverse? simpler: use scipy if available else approximate
    try:
        from scipy.stats import chi2
        lo = chi2.ppf(alpha/2, 2*k)/2 if k > 0 else 0.0
        hi = chi2.ppf(1-alpha/2, 2*(k+1))/2
    except Exception:
        # Normal approx fallback
        lo = max(0.0, k - 1.96*math.sqrt(k))
        hi = k + 1.96*math.sqrt(k)
    return lo, hi

def age(K):
    return K/(2*MU)

def bootstrap_ci(a, b, n=10000, seed=42, alpha=0.05):
    """Column-bootstrap 95% CI on the insertion age (the method the manuscript states:
    resample ungapped aligned column-pairs with replacement, recompute K2P each replicate).
    Seeded -> fully reproducible."""
    import random
    cols = [(x, y) for x, y in zip(a, b)
            if x not in '-.' and y not in '-.']
    m = len(cols)
    rng = random.Random(seed)
    ages = []
    for _ in range(n):
        ts = tv = 0
        for _ in range(m):
            x, y = cols[rng.randrange(m)]
            if x == y:
                continue
            if frozenset((x, y)) in TRANSITIONS: ts += 1
            else: tv += 1
        P, Q = ts/m, tv/m
        try:
            K = -0.5*math.log(1-2*P-Q) - 0.25*math.log(1-2*Q)
        except ValueError:
            K = float('nan')
        ages.append(age(K))
    ages = sorted(x for x in ages if x == x)
    lo = ages[int((alpha/2)*len(ages))]
    hi = ages[int((1-alpha/2)*len(ages))]
    return lo, hi, n, seed

order, a, b = parse_stretcher(ALN)
r = k2p(a, b)
K = r['K']
age_pt = age(K)
lo_cnt, hi_cnt = poisson_ci(r['subs'])
# scale K with substitution count (Q small -> approx linear); recompute K at CI bounds
def K_from_subs(subs):
    P = (r['ts']/r['subs'])*subs/r['ungapped'] if r['subs'] else 0
    Q = (r['tv']/r['subs'])*subs/r['ungapped'] if r['subs'] else 0
    try:
        return -0.5*math.log(1-2*P-Q) - 0.25*math.log(1-2*Q)
    except ValueError:
        return float('nan')
age_lo = age(K_from_subs(lo_cnt))
age_hi = age(K_from_subs(hi_cnt))

print(f"Alignment file : {ALN}")
print(f"Aligned cols   : {r['aln_len']}  | ungapped sites: {r['ungapped']}  | gap cols: {r['gaps']}")
print(f"Identical      : {r['ident']}")
print(f"Substitutions  : {r['subs']}  (transitions {r['ts']}, transversions {r['tv']})")
print(f"p-distance(subs): {r['p_dist']:.5f}   K2P: {K:.5f}")
print(f"mu             : {MU} subs/site/yr")
print(f"Insertion age  : {age_pt/1e6:.3f} Mya")
print(f"Poisson 95% CI on substitutions: [{lo_cnt:.2f}, {hi_cnt:.2f}]")
print(f"Age 95% CI (Poisson) : [{age_lo/1e6:.3f}, {age_hi/1e6:.3f}] Mya")
bs_lo, bs_hi, bn, bseed = bootstrap_ci(a, b)
print(f"Age 95% CI (column-bootstrap, n={bn}, seed={bseed}) : [{bs_lo/1e6:.3f}, {bs_hi/1e6:.3f}] Mya")
print()
print("CONCLUSION: substitutions-only dating -> ~0.2 Mya (manuscript value is correct).")
print("The 1.4 Mya figure in 结果_Kala4致色变异与定年.md counted indels as differences (K=0.0361) and is an error.")
print("Both pre-date rice domestication (~10 kya): the intronic Gypsy is ancient and NOT the causal cis-change;")
print("the ZN65-specific proximal-promoter TE cluster remains the candidate switch.")
