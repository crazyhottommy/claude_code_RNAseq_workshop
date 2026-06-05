#!/usr/bin/env bash
#
# 03_salmon_quant.sh
# ------------------
# Quantify each single-end FASTQ against the salmon index, producing one
# quant.sf per sample under data/salmon/<sample>/.
#
# Usage:
#   bash scripts/03_salmon_quant.sh
#
# Requires the `rnaseq` conda env (salmon) and a built index (script 02).
set -euo pipefail

INDEX="${INDEX:-reference/gencode.v45_human_index}"
FASTQDIR="${FASTQDIR:-fastq}"
OUTROOT="${OUTROOT:-data/salmon}"
THREADS="${THREADS:-4}"

SAMPLES=(Normoxia_sgCTRL_1 Normoxia_sgCTRL_2 Hypoxia_sgCTRL_1 Hypoxia_sgCTRL_2)

for name in "${SAMPLES[@]}"; do
  echo ">>> quantifying $name"
  # -l A      : auto-detect library type
  # -r        : single-end reads (paired-end would use -1 / -2)
  # --validateMappings : selective-alignment for higher accuracy
  # --gcBias / --seqBias : correct common technical biases
  salmon quant \
    -i "$INDEX" \
    -l A \
    -r "$FASTQDIR/${name}.fastq.gz" \
    -o "$OUTROOT/${name}" \
    --validateMappings \
    --gcBias \
    --seqBias \
    -p "$THREADS"
done

echo "Quantification done. quant.sf files under $OUTROOT/"
