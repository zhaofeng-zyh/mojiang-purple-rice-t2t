# Mojiang purple rice (ZN65) T2T genome — analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21236668.svg)](https://doi.org/10.5281/zenodo.21236668)


Analysis and figure-generation code accompanying the manuscript:

> **A gap-free telomere-to-telomere genome of the purple-pericarp rice landrace Mojiang ZN65 resolves transposon-associated structural variation at the pericarp-pigmentation loci.**
> Feng Zhao, Juan Zhao, Fang Zhao, Suhong Bai, Yahan Wu, Rengguo Zhu. Pu'er University. (submitted, 2026)

This repository contains the **analysis pipeline** and the **small derived tables** it produces/consumes. It does **not** contain the genome assembly, annotation, or raw sequencing reads (see *Data availability*).

## Repository layout

```
pipeline/
├── steps/                # one script per analysis stage (as-run on the authors' server)
│   ├── a3_longread.sh / a3_hifi_analyze.sh     # long-read validation of pigmentation-locus SVs
│   ├── a4_align.sh / a4_alleles.sh             # functional-allele status of the four loci
│   ├── a9_sv.sh / a9_compare.sh / a9_assemblytics.sh   # genome-wide SV (SyRI/minimap2) + orthogonal calls
│   ├── genome_landscape.sh                     # gene/TE/GC windows for the Circos landscape
│   ├── a11_pathway_rbh.sh                       # flavonoid-pathway inventory (KEGG/RBH)
│   ├── a12_*.sh / a12_*.py                      # 3K-RGP subspecies placement, gene trees
│   ├── a14_wd40.sh                             # OsTTG1 (WD40) identification
│   ├── a16_te.sh                              # transposable-element family attribution
│   ├── of_*.sh / batch2*_qc.sh                # orthology, assembly QC
│   ├── centromere_rDNA.py                     # CentO centromeres + 45S/5S rDNA arrays  (Table S5)
│   ├── pathway_copynumber.py                  # ZN65-vs-Nipponbare pathway copy number  (Table S4, Fig 3A)
│   ├── expression.R                           # TPM + DESeq2 over vegetative RNA-seq     (Table S12, Fig S3)
│   ├── ltr_dating.py                          # K2P dating of the intronic Gypsy LTRs    (~0.2 Myr)
│   ├── kala1_helitron_dup.py                  # Helitron duplication at Kala1/OsDFR
│   └── build_supp_tables.py                   # assembles the 13-sheet Supplementary_Tables.xlsx
└── environment/         # R and Python dependency lists
data/                    # small derived tables consumed/produced by the scripts (no raw data)
```

## Dependencies

- **R ≥ 4.4** — `circlize`, `ggplot2`, `patchwork`, `ggtree`, `treeio`, `ComplexHeatmap`, `pheatmap`, `DESeq2`, `data.table` (`pipeline/environment/R_packages.txt`).
- **Python ≥ 3.10** — `matplotlib`, `numpy`, `pandas`, `openpyxl`, `python-docx` (`pipeline/environment/python_requirements.txt`).
- **External tools (upstream, versions per the paper Methods):** hifiasm 0.25.0, NextPolish 1.4.1, Juicer 1.6 + 3D-DNA, quarTeT 1.2.5, Merqury 1.3, LTR_retriever, BUSCO 5.8.0, RepeatModeler2/RepeatMasker, minimap2, SyRI 1.7.1, plotsr, MUMmer4, MAFFT, IQ-TREE 2, OrthoFinder 2.5, HISAT2 + featureCounts, PLINK 1.9/2.0.

## Reproducing the analyses

The `pipeline/steps/` scripts were written to run on the authors' compute environment and contain **machine-specific absolute paths**. They are released as an accurate record of the analyses, not as a turnkey workflow — adjust paths to your local copies of the assembly, annotation and reads. The `data/` folder provides the small input and derived-output tables used by these analyses (genome windows, SyRI summary, LTR alignment, pathway list, RNA-seq count matrix, and the result tables).

### Adapting the paths

`adapt_paths.sh` rewrites those absolute paths to your own layout in one step:

```bash
bash adapt_paths.sh --check            # report affected files, change nothing
bash adapt_paths.sh /data/my_zn65      # rewrite in place (timestamped backup is made first)
```

It currently affects **89 files / 156 lines**. After rewriting it runs `bash -n` and `python3 -m py_compile`
over every modified script and refuses to finish if any of them stops parsing; the printed backup path
restores the original state.

Your target directory has to mirror the layout the scripts expect, for example
`<root>/01_基因组_组装与注释_Genome/01.assembly/ZN65.T2T.fa` and `<root>/07_分析_Os02g基因鉴定/…`.
The sequence data themselves are not in this repository — see **Data availability** below.

> **Why literal substitution rather than a `${ROOT}` variable.** Ten of these scripts embed Python
> through *quoted* here-documents (`<<'PY'`), six of which also contain absolute paths. Inside a quoted
> here-document the shell performs no parameter expansion, so a `${ROOT}`-style rewrite would hand the
> literal text `${ROOT}` to Python and break those steps silently. Literal-to-literal replacement is
> immune to that failure mode, which is why it is done this way.

⚠️ Some scripts under `pipeline/steps/` are **superseded** and fail closed on purpose (they reproduce
figures under an interpretation that has since been corrected). They print an explanation and exit
non-zero unless `ALLOW_SUPERSEDED=1` is set; output produced that way must not be used as current results.

## Data availability

- **BioProject:** PRJCA068452 (CNCB / NGDC).
- **Raw reads (GSA):** submission `subCRA075071` — PacBio HiFi, ONT ultra-long, DNBSEQ WGS, Hi-C and RNA-seq runs `CRR3301760`–`CRR3301769`.
- **Genome assembly + annotation (GWH):** batch `Batch0093257` (accession `GWH…`, released on publication).
- **Reference genomes:** Nipponbare IRGSP-1.0; MH63RS3/ZS97RS3; Cempo Ireng GCA_055776245.1; N22 GCA_001952365.2; *O. rufipogon* W1943 GCA_000817225.1; 3,000 Rice Genomes Project.

`data/deposit_md5_manifest.tsv` lists the deposited raw files with MD5 checksums.

## Citation

If you use this code, please cite the manuscript above and the archived code:

> Zhao, F., Zhao, J., Zhao, F., Bai, S., Wu, Y., Zhu, R. (2026). Mojiang purple rice (ZN65) T2T genome — analysis code. Zenodo. https://doi.org/10.5281/zenodo.21236668

## License

Released under the MIT License (see [LICENSE](LICENSE)).
