#!/usr/bin/env bash
#
# 01_download_fastq.sh
# --------------------
# Download the 4 single-end FASTQ files for GEO series GSE197576.
#
#   Normoxia sgCTRL : SRX14311105, SRX14311106
#   Hypoxia  sgCTRL : SRX14311111, SRX14311112
#
# Two modes:
#   FULL   (default) - fastq-dl pulls the complete FASTQ from ENA. Faithful to
#                      the blog post, but several GB per sample.
#   SUBSET           - sra-tools fastq-dump streams only the first N reads.
#                      Fast and small; this is how the quant.sf committed to
#                      this repo were produced (NREADS=3,000,000).
#
# Usage:
#   bash scripts/01_download_fastq.sh            # FULL
#   MODE=SUBSET bash scripts/01_download_fastq.sh
#
# Requires the `rnaseq` conda env (fastq-dl, sra-tools). See environment.yml.
set -euo pipefail

MODE="${MODE:-FULL}"
NREADS="${NREADS:-3000000}"
OUTDIR="${OUTDIR:-fastq}"
mkdir -p "$OUTDIR"

# SRX experiment accession -> human-readable sample name used everywhere downstream.
declare -a SAMPLES=(
  "SRX14311105:Normoxia_sgCTRL_1"
  "SRX14311106:Normoxia_sgCTRL_2"
  "SRX14311111:Hypoxia_sgCTRL_1"
  "SRX14311112:Hypoxia_sgCTRL_2"
)

for entry in "${SAMPLES[@]}"; do
  srx="${entry%%:*}"
  name="${entry##*:}"
  echo ">>> $name ($srx)  mode=$MODE"

  if [[ "$MODE" == "SUBSET" ]]; then
    # fastq-dump resolves the SRX -> SRR and streams only the first NREADS spots.
    fastq-dump -X "$NREADS" --gzip -O "$OUTDIR" "$srx"
    # fastq-dump names the file after the SRR run; rename to our sample name.
    newest=$(ls -t "$OUTDIR"/*.fastq.gz | head -1)
    mv "$newest" "$OUTDIR/${name}.fastq.gz"
  else
    # fastq-dl groups all runs of an experiment into one FASTQ.
    fastq-dl --accession "$srx" --group-by-experiment --outdir "$OUTDIR"
    newest=$(ls -t "$OUTDIR"/*.fastq.gz | head -1)
    mv "$newest" "$OUTDIR/${name}.fastq.gz"
  fi
done

echo "Done. FASTQ in $OUTDIR/"
