#!/usr/bin/env bash
#
# 02_build_salmon_index.sh
# ------------------------
# Download the GENCODE v45 human transcriptome and build a salmon index.
#
# This is the compute/RAM-heavy step (a transcriptome index needs a few GB of
# RAM and ~10-20 min). It is NOT run during the live workshop -- the committed
# quant.sf already exist. Run it only if you want to re-quantify from FASTQ.
#
# Usage:
#   bash scripts/02_build_salmon_index.sh
#
# Requires the `rnaseq` conda env (salmon). See environment.yml.
set -euo pipefail

REFDIR="${REFDIR:-reference}"
INDEX="${INDEX:-reference/gencode.v45_human_index}"
THREADS="${THREADS:-4}"
TX_FA="$REFDIR/gencode.v45.transcripts.fa.gz"
GTF="$REFDIR/gencode.v45.annotation.gtf.gz"
BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45"

mkdir -p "$REFDIR"

# Transcript FASTA (used to build the index).
[[ -f "$TX_FA" ]] || wget -O "$TX_FA" "$BASE/gencode.v45.transcripts.fa.gz"
# GTF annotation (used later by make_tx2gene.R to map transcript -> gene).
[[ -f "$GTF" ]]   || wget -O "$GTF"   "$BASE/gencode.v45.annotation.gtf.gz"

# --gencode strips the trailing "|"-delimited metadata from GENCODE FASTA headers
# so the index uses clean ENST IDs. -k 31 is the standard k-mer size for >=75bp reads.
salmon index \
  -t "$TX_FA" \
  -i "$INDEX" \
  -k 31 \
  --gencode \
  -p "$THREADS"

echo "Index built at $INDEX"
