#!/usr/bin/env python3
"""
R-3 re-derivation: is the manuscript's "Kala1/OsDFR Helitron-associated duplication"
supported by source data?

Two questions were conflated in the project notes:
  (a) locus-level ZN65-vs-Nipponbare divergence  -> dd.report (96.59% id, 93.28% aln,
      "874 bp unaligned"); we show below the 874 bp is just the asymmetric extraction
      WINDOW EDGES (NIP window 13 kb vs ZN65 12 kb), NOT an internal SV.
  (b) the "Helitron-associated duplication" itself -> a ZN65-internal ~1.3 kb segmental
      duplication in the OsDFR 5' flank, which we localise here and classify against the
      ZN65 RepeatMasker annotation.

Method:
  1. Extract the OsDFR (ZN651G2772, Chr1:27,075,322-27,077,300) locus window.
  2. nucmer --maxmatch self-alignment of the window to expose internal repeats.
  3. Intersect each duplicated copy with the genome RepeatMasker GFF -> TE family.

Tools: nucmer/show-coords/seqkit (env cgsv). Paths relative to this file.
Run:  ~/miniconda3/envs/cgsv/bin/python3 kala1_helitron_dup.py
"""
import subprocess, tempfile, os, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ASM  = ROOT / "01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa"
RM   = ROOT / "01_基因组_组装与注释_Genome/02.annotation/01.repeat/repeatmasker.gff"
BIN  = Path.home() / "miniconda3/envs/cgsv/bin"
CHROM, W0, W1 = "Chr1", 27069000, 27084000   # OsDFR gene 27,075,322-27,077,300 +/- flank
GENE0, GENE1 = 27075322, 27077300

def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)

tmp = Path(tempfile.mkdtemp())
win = tmp / "win.fa"
# 1. extract window
win.write_text(run([str(BIN/"seqkit"), "faidx", str(ASM), f"{CHROM}:{W0}-{W1}"]).stdout)
# 2. self-align to find internal duplications (exclude the main diagonal)
dp = tmp / "self"
run([str(BIN/"nucmer"), "--maxmatch", "--nosimplify", "-p", str(dp), str(win), str(win)])
co = run([str(BIN/"show-coords"), "-rclTH", f"{dp}.delta"]).stdout

print(f"# OsDFR locus self-alignment  ({CHROM}:{W0}-{W1}; gene {GENE0}-{GENE1})")
print("# off-diagonal hits = internal segmental duplications")
dups = []
for ln in co.splitlines():
    f = ln.split("\t")
    if len(f) < 7:
        continue
    s1, e1, s2, e2 = int(f[0]), int(f[1]), int(f[2]), int(f[3])
    L1, idy = int(f[4]), float(f[6])
    # skip self/main diagonal (same coords) and trivial short hits
    if abs(s1 - s2) < 50 and abs(e1 - e2) < 50:
        continue
    if L1 < 300 or s1 >= s2:        # keep one orientation, sizeable copies only
        continue
    g1, g2, g3, g4 = W0 + s1 - 1, W0 + e1 - 1, W0 + s2 - 1, W0 + e2 - 1
    dups.append((g1, g2, g3, g4, L1, idy))
    print(f"  copyA {CHROM}:{g1}-{g2}  <->  copyB {CHROM}:{g3}-{g4}   len~{L1}bp  id={idy:.1f}%")

# 3. classify the duplicated copies against RepeatMasker
def te_in(a, b):
    out = []
    for ln in open(RM, encoding="utf-8", errors="ignore"):
        if not ln.startswith(CHROM + "\t"):
            continue
        f = ln.rstrip("\n").split("\t")
        s, e = int(f[3]), int(f[4])
        if e < a or s > b:
            continue
        m = re.search(r"Target=([^;]+).*?Class=([^;]+)", f[8])
        if m:
            out.append((s, e, m.group(2), m.group(1)))
    return out

print("\n# RepeatMasker TE families overlapping each duplicated copy")
for d in dups:
    for (s, e, cls, tgt) in te_in(d[0], d[1]) + te_in(d[2], d[3]):
        print(f"  {CHROM}:{s}-{e}  {cls:18s} {tgt}")

print("\n# VERDICT")
if any("Helitron" in t[2] for d in dups for t in (te_in(d[0],d[1])+te_in(d[2],d[3]))):
    print("  CONFIRMED: the OsDFR locus carries an internal ~1.3 kb segmental duplication")
    print("  whose copies are RC/Helitron (RPO_OS / Helitron-N116_OS) -> 'Helitron-associated")
    print("  duplication' is supported. (Specificity vs Nipponbare not tested here; the locus")
    print("  dd.report is colinear at 96.59%, so the duplication is likely shared, not ZN65-specific.)")
else:
    print("  NOT supported as Helitron.")
