# install.R
# ---------
# Install the R / Bioconductor packages used in the downstream analysis.
# Run once:  Rscript install.R

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cran <- c("dplyr", "readr", "tibble", "ggplot2", "pheatmap", "ggrepel")
bioc <- c("tximport", "DESeq2", "rtracklayer", "clusterProfiler",
          "org.Hs.eg.db", "AnnotationDbi")

to_install <- setdiff(c(cran, bioc), rownames(installed.packages()))
if (length(to_install)) {
  BiocManager::install(to_install, update = FALSE, ask = FALSE)
} else {
  message("All packages already installed.")
}
