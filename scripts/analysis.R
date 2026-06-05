#!/usr/bin/env Rscript
# analysis.R
# ----------
# Complete reference downstream RNAseq analysis for the workshop:
#   salmon quant.sf  ->  tximport  ->  DESeq2  ->  results + figures
#
# This is the "instructor copy" / safety net. During the live workshop the same
# steps are produced interactively by prompting Claude Code (see the website).
#
# Usage:  Rscript scripts/analysis.R
# Outputs: results/de_hypoxia_vs_normoxia.csv and results/figures/*.png

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
save_plot <- function(p, file, w = 6, h = 5) {
  ggsave(file.path("results/figures", file), p, width = w, height = h, dpi = 150)
}

# ---- 1. Sample sheet -------------------------------------------------------
samples <- read_csv("data/samples.csv", show_col_types = FALSE)
samples$condition <- factor(samples$condition, levels = c("normoxia", "hypoxia"))
files <- file.path("data", "salmon", samples$sample, "quant.sf")
names(files) <- samples$sample
stopifnot(all(file.exists(files)))

# ---- 2. tximport: transcript -> gene counts --------------------------------
tx2gene <- read_csv("data/tx2gene.csv", show_col_types = FALSE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene,
                ignoreTxVersion = FALSE)

# ---- 3. DESeq2 dataset -----------------------------------------------------
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "normoxia")

# ---- 4. Filter low-count genes --------------------------------------------
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]
message("Genes after filtering: ", nrow(dds))

# ---- 5. Differential expression -------------------------------------------
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "hypoxia", "normoxia"))

# ---- 6. Annotate with gene symbols ----------------------------------------
gene_name_map <- read_csv("data/gene_name_map.csv", show_col_types = FALSE)
res_df <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_name_map, by = "gene_id") %>%
  arrange(padj)

dir.create("results", showWarnings = FALSE)
write_csv(res_df, "results/de_hypoxia_vs_normoxia.csv")
# Saved so the website's visualization/enrichment pages can load a ready DESeq2
# object without recomputing (and independent of page render order).
saveRDS(dds, "results/dds.rds")

sig <- res_df %>% filter(!is.na(padj) & padj < 0.05)
message("Significant genes (padj < 0.05): ", nrow(sig))
message("Top hits: ", paste(head(na.omit(sig$gene_name), 10), collapse = ", "))

# ---- 7. Figures ------------------------------------------------------------
# p-value histogram
ph <- ggplot(res_df, aes(pvalue)) +
  geom_histogram(boundary = 0, bins = 30, fill = "#3a6ea5", color = "white") +
  labs(title = "p-value distribution", x = "p-value", y = "count") +
  theme_bw()
save_plot(ph, "pvalue_histogram.png")

# PCA on variance-stabilized counts
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"))
pca_p <- ggplot(pca, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  labs(x = paste0("PC1: ", pv[1], "% var"),
       y = paste0("PC2: ", pv[2], "% var"), title = "PCA") +
  theme_bw()
save_plot(pca_p, "pca.png")

# MA plot
png("results/figures/ma_plot.png", width = 900, height = 750, res = 150)
plotMA(res, ylim = c(-5, 5), main = "MA plot: hypoxia vs normoxia")
dev.off()

# Volcano
volc <- res_df %>%
  mutate(sig = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
volc_p <- ggplot(volc, aes(log2FoldChange, -log10(pvalue), color = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#c0392b"),
                     guide = "none") +
  labs(title = "Volcano: hypoxia vs normoxia",
       x = "log2 fold change", y = "-log10 p-value") +
  theme_bw()
save_plot(volc_p, "volcano.png")

# Heatmap of top 30 DE genes (vst, z-scored by row)
top_ids <- head(sig$gene_id, 30)
if (length(top_ids) >= 2) {
  mat <- assay(vsd)[top_ids, , drop = FALSE]
  rownames(mat) <- gene_name_map$gene_name[match(top_ids, gene_name_map$gene_id)]
  ann <- data.frame(condition = samples$condition)
  rownames(ann) <- samples$sample
  pheatmap(mat, scale = "row", annotation_col = ann,
           show_rownames = TRUE, fontsize_row = 7,
           filename = "results/figures/heatmap_top30.png",
           width = 7, height = 7)
}

message("Done. See results/de_hypoxia_vs_normoxia.csv and results/figures/")
