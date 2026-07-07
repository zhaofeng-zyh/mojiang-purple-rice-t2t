import re,sys
BM="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/V2/bm.py"
s=open(BM,encoding="utf-8").read()
caps={
"Fig1_genome_landscape.png":
 "Telomere-to-telomere genome landscape of purple rice ZN65 (indica) across twelve gap-free chromosomes "
 "(395.1 Mb). Rings, outer to inner: chromosome ideogram (10-Mb tick labels; orange dots, resolved telomeres "
 "[24/24]; stars with leaders, the four pericarp-pigmentation loci OsB2/Kala4, OsDFR/Kala1, OsC1/Kala3, Rc; "
 "triangles, approximate centromere taken as the maximum-TE 500-kb window); gene density (genes per 500 kb); "
 "transposable-element coverage (%); GC content (line, 41-47%). Heatmap intensities are globally 2-98-percentile "
 "scaled per track (comparable across chromosomes); densities and coverage were computed in 500-kb windows with "
 "bedtools. 45S/5S rDNA arrays were not separately annotated in this version.",
"Fig2_genomewide_SV.png":
 "Genome-wide structural divergence and transposon-driven expansion of ZN65 versus Nipponbare. (a) Whole-genome "
 "synteny (minimap2 asm5 + SyRI, visualised with plotsr): grey, syntenic; orange, inversions (note the large "
 "chromosome-6 inversion); green, translocations; blue, duplications. (b) Variant spectrum (log scale) split into "
 "sequence-level (SNP, insertion, deletion) and structural (inversion, translocation, duplication, highly diverged) "
 "classes; the 449 translocations quoted in text are the SyRI aggregate of TRANS (198) and inverted translocations "
 "(251). (c) Transposable-element families and non-TE unique sequence within the 73.3 Mb of ZN65-specific (NOTAL+INS) "
 "sequence (RepeatMasker). (d) Size distribution of structural variants (INV/TRANS/DUP). Robustness to minimap2 "
 "stringency (asm5/10/20) and an orthogonal Assemblytics call are reported in Methods.",
"Fig3_pathway.png":
 "Flavonoid/anthocyanin biosynthetic pathway in ZN65 versus Nipponbare. Each enzyme box shows copy number as "
 "ZN65 / Nipponbare, assigned by KEGG orthology and counted as reciprocal-best-hit orthologues (may include "
 "paralogues/pseudogenes; high-confidence curation and an OrthoFinder orthogroup baseline are pending). Intermediate "
 "metabolites are labelled. The pathway is complete and copy-number-conserved between the two genomes, so the "
 "pigmentation difference maps to the MBW regulators (shown as a model), not the enzymes. The anthocyanin end-products "
 "(cyanidin- and peonidin-3-O-glucoside, C3G/P3G) are not yet quantified and are flagged 'to be quantified "
 "(HPLC-DAD, G1)'.",
"Fig3_loci.png":
 "Gene-model and transposable-element architecture of the four pericarp-pigmentation loci in ZN65, showing where "
 "structural variation falls relative to gene structure. For each locus the gene model (blue exons, intron line, "
 "strand arrow) is shown above its RepeatMasker TE annotation (coloured by class). Kala4/OsB2 carries a ~12.5-kb "
 "Gypsy within a 14.4-kb intron plus a ZN65-specific promoter insertion; Rc carries Stowaway/MITEs in its promoter "
 "and introns; Kala1/OsDFR has a Helitron cluster and Stowaway flanking the promoter; OsC1/Kala3 is colinear and "
 "conserved (99.8% identity), with TEs only in flanks. Per-panel scale differs (1-kb bar in each). Models from the "
 "ZN65 annotation (GFF) and RepeatMasker.",
"Fig4_kala4_architecture.png":
 "Retrotransposon architecture and lineage-specific origin of the Kala4/OsB2 pigmentation allele in ZN65 "
 "(chromosome 4, minus strand). (a) OsB2 gene model with exon-intron structure (auto-read from the annotation): a "
 "complete Gypsy LTR retrotransposon (RETRO2B, ~12.5 kb) sits within the largest (~14.4-kb) intron and a 6.3-kb "
 "LINE-1 within the gene, with a nested promoter TE cluster; the ZN65-specific promoter insertion spans ~5.9 kb. The "
 "intronic Gypsy LTRs date to approximately 0.2 Myr (K2P, 95% CI 0.04-0.39; mu = 1.3e-8). (b) The ZN65-specific "
 "insertion, scored by local coverage within the orthologous OsB2 locus (anchored on the 34-kb conserved region), is "
 "present only in ZN65 and absent from indica (MH63, ZS97), aus (N22), wild (O. rufipogon) and the black-pericarp "
 "landrace Cempo Ireng. Caveat: one genome per accession/species; intraspecific polymorphism was not assessed.",
}
def setcap(s,fname,newcap):
    pat=re.compile(r'(FIGURE\("'+re.escape(fname)+r'",\s*"[^"]*",\s*)(?:"(?:[^"\\]|\\.)*"\s*)+(\))',re.S)
    n=len(pat.findall(s))
    if n!=1: print(f"ERR {fname} matches={n}"); sys.exit(1)
    lit='"'+newcap.replace('"',"'")+'"'
    return pat.sub(lambda m:m.group(1)+lit+m.group(2),s,count=1)
for fn,cap in caps.items():
    s=setcap(s,fn,cap); print("caption set:",fn)
open(BM,"w",encoding="utf-8").write(s)
print("all 4 captions updated")
