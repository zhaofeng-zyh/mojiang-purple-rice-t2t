suppressMessages({library(data.table)})
BASE <- "/mnt/data2/墨江紫米研究"
OUT  <- file.path(BASE,"12_论文Paper1_Manuscript/_SUBMISSION_PlantComm/analysis")
# --- counts + gene lengths (from featureCounts) ---
fc <- fread(file.path(BASE,"14_转录组_重比对ZN65T2T/counts/ZN65_gene_counts.txt"), skip=1)
setnames(fc, 1, "Geneid")
samples <- c("SD1","SD2","SD3","TG1","TG3")
# featureCounts columns end with the bam paths; map them
bamcols <- grep("\\.sorted\\.bam$", names(fc), value=TRUE)
stopifnot(length(bamcols)==5)
cnt <- as.matrix(fc[, ..bamcols]); rownames(cnt) <- fc$Geneid; colnames(cnt) <- samples
len <- fc$Length; names(len) <- fc$Geneid
# --- TPM ---
rpk <- cnt / (len/1000)
tpm <- t(t(rpk) / (colSums(rpk)/1e6))
# --- key genes ---
key <- data.table(
  gene = c("ZN656G0716","ZN654G2687","ZN657G0823","ZN652G3195","ZN651G2772"),
  symbol = c("OsC1/Kala3 (R2R3-MYB)","OsB2/Kala4 (bHLH)","Rc (bHLH)","OsTTG1 (WD40)","OsDFR/Kala1 (DFR)"),
  role = c("MBW-MYB","MBW-bHLH(anthocyanin)","bHLH(proanthocyanidin)","MBW-WD40","structural"))
kt <- cbind(key, round(tpm[key$gene, , drop=FALSE],2))
cat("=== Key pigmentation/MBW gene TPM (vegetative tissue: SD=seedling n=3, TG=tillering n=2; TG2 absent) ===\n")
print(kt)
# means per condition
kt[, SD_mean := round(rowMeans(tpm[key$gene, c("SD1","SD2","SD3")]),2)]
kt[, TG_mean := round(rowMeans(tpm[key$gene, c("TG1","TG3")]),2)]
fwrite(kt, file.path(OUT,"expression_key_TPM.tsv"), sep="\t")

# --- pathway genes TPM summary ---
gp <- readLines(file.path(BASE,"07_分析_Os02g基因鉴定/18_通路拷贝数_A11/genes.txt"))
gp <- gp[gp %in% rownames(tpm)]
pw_tpm <- tpm[gp,,drop=FALSE]
pw_mean <- rowMeans(pw_tpm)
cat(sprintf("\n=== Pathway genes (n=%d of 100 found in count matrix): median TPM=%.2f; %% with mean TPM>=1: %.0f%% ===\n",
    length(gp), median(pw_mean), 100*mean(pw_mean>=1)))
fwrite(data.table(gene=gp, round(pw_tpm,2), mean_TPM=round(pw_mean,2)), file.path(OUT,"expression_pathway_TPM.tsv"), sep="\t")

# --- DESeq2 stabilized view (SD vs TG) - flagged low power ---
suppressMessages(library(DESeq2))
coldata <- data.frame(row.names=samples, tissue=factor(c("seedling","seedling","seedling","tillering","tillering")))
dds <- DESeqDataSetFromMatrix(round(cnt), coldata, ~tissue)
dds <- dds[rowSums(counts(dds))>=10,]
dds <- DESeq(dds, quiet=TRUE)
res <- as.data.frame(results(dds, contrast=c("tissue","tillering","seedling")))
res$gene <- rownames(res)
cat(sprintf("\n=== DESeq2 SD vs TG (CAVEAT: n=3 vs n=2, TG2 missing; exploratory only) ===\nGenes tested: %d; padj<0.05: %d\n",
    nrow(res), sum(res$padj<0.05, na.rm=TRUE)))
# key genes DESeq stats
print(round(res[key$gene, c("baseMean","log2FoldChange","pvalue","padj")],3))
fwrite(as.data.table(res)[order(padj)], file.path(OUT,"expression_DESeq2_SDvsTG.tsv"), sep="\t")
saveRDS(list(tpm=tpm, key=kt), file.path(OUT,"expression_objs.rds"))
cat("\nWrote expression_key_TPM.tsv, expression_pathway_TPM.tsv, expression_DESeq2_SDvsTG.tsv\n")
