# CLAUDE.md — project context for Claude Code

This repo is a hands-on workshop: **analyze bulk RNAseq data end-to-end using Claude Code**.
When an attendee opens Claude Code here, use this context to help them.

## The dataset
- GEO series **GSE197576**, human bulk RNAseq, single-end reads.
- 4 samples, 2 conditions (2 replicates each):
  - `normoxia`: `Normoxia_sgCTRL_1` (SRX14311105), `Normoxia_sgCTRL_2` (SRX14311106)
  - `hypoxia` : `Hypoxia_sgCTRL_1` (SRX14311111),  `Hypoxia_sgCTRL_2` (SRX14311112)
- Biological expectation: hypoxia induces HIF targets (e.g. **VEGFA, CA9, BNIP3, SLC2A1/GLUT1, PGK1**).

## Where things live
- `data/salmon/<sample>/quant.sf` — salmon transcript quantifications (the input to load).
- `data/samples.csv` — sample sheet (`sample,condition,srx,replicate`).
- `data/tx2gene.csv` — transcript_id → gene_id (for tximport gene-level summarization).
- `data/gene_name_map.csv` — gene_id → gene_name (for annotating results).
- `scripts/` — preprocessing shell scripts + `analysis.R` (full reference solution).
- `*.qmd` — the Quarto website pages (the workshop guide).

## The analysis we want (downstream)
1. Load salmon output with **tximport** (`type = "salmon"`, `tx2gene`).
2. Build a **DESeq2** dataset: `DESeqDataSetFromTximport(txi, colData, ~ condition)`,
   `relevel(condition, ref = "normoxia")`.
3. Filter low-count genes (`rowSums(counts) > 10`), run `DESeq()`.
4. `results(dds, contrast = c("condition", "hypoxia", "normoxia"))`, annotate symbols via `gene_name_map`.
5. Visualize: p-value histogram, PCA (vst), MA plot, volcano, heatmap of top genes.
6. Functional enrichment with **clusterProfiler::enrichGO** on significant genes.

## Conventions
- R + Bioconductor. Preferred packages: tximport, DESeq2, tidyverse (dplyr/readr/ggplot2), pheatmap, clusterProfiler, org.Hs.eg.db.
- Gene IDs are **versioned GENCODE** ENSG/ENST (e.g. `ENSG00000....17`); keep versions consistent between tx2gene and quant.sf.
- Save figures under `results/figures/`, tables under `results/`.
- The reference end-to-end script is `scripts/analysis.R` — match its logic if asked to reproduce results.

## Not run live
Preprocessing (`scripts/01–03`, salmon index build, FASTQ download) is documented but not run during the 1-hour session; the committed `quant.sf` are the starting point.
