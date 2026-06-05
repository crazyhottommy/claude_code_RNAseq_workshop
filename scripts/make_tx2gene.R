#!/usr/bin/env Rscript
# make_tx2gene.R
# --------------
# Build the transcript->gene map (tx2gene) and a gene_id->gene_name map from the
# GENCODE v45 GTF. tximport uses tx2gene to aggregate transcript counts to genes;
# gene_name_map is used downstream to label results with HGNC symbols.
#
# Usage:  Rscript scripts/make_tx2gene.R [path/to/gencode.v45.annotation.gtf.gz]
#
# Writes: data/tx2gene.csv, data/gene_name_map.csv

suppressPackageStartupMessages({
  library(rtracklayer)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
gtf_path <- if (length(args) >= 1) args[1] else "reference/gencode.v45.annotation.gtf.gz"

if (!file.exists(gtf_path)) {
  stop("GTF not found at '", gtf_path, "'. Run scripts/02_build_salmon_index.sh first ",
       "(it downloads gencode.v45.annotation.gtf.gz), or pass the path as an argument.")
}

message("Importing GTF: ", gtf_path)
gtf <- rtracklayer::import(gtf_path)
gtf_df <- as.data.frame(gtf)

# Transcript-level rows carry both transcript_id and gene_id.
tx <- gtf_df %>%
  filter(type == "transcript") %>%
  select(transcript_id, gene_id, gene_name) %>%
  distinct()

dir.create("data", showWarnings = FALSE)

tx %>%
  select(transcript_id, gene_id) %>%
  write_csv("data/tx2gene.csv")

tx %>%
  select(gene_id, gene_name) %>%
  distinct() %>%
  write_csv("data/gene_name_map.csv")

message("Wrote data/tx2gene.csv (", nrow(tx), " transcripts) and data/gene_name_map.csv")
