#!/bin/bash
# Quantifies every library, reverse-stranded. Skips anything already done.
set -u
K="$HOME/Downloads/kallisto/kallisto"; IDX="$HOME/Downloads/index.idx"
BASE="$HOME/Library/CloudStorage/OneDrive-pnw.edu/illumina basespace content/BaseSpace"
OUT="$HOME/postn_rnaseq/quant"; mkdir -p "$OUT"
find "$BASE" -name "*_R1_001.fastq.gz" -print0 | while IFS= read -r -d '' R1; do
  R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
  NAME=$(basename "$R1" | sed -E 's/_S[0-9]+_L[0-9]+_R1_001\.fastq\.gz$//')
  [ -s "$OUT/$NAME/abundance.tsv" ] && { echo "$NAME done"; continue; }
  "$K" quant -i "$IDX" -o "$OUT/$NAME" -t 4 --rf-stranded "$R1" "$R2"
done
for d in "$OUT"/*/; do
  n=$(basename "$d")
  p=$(grep -o '"p_pseudoaligned": [0-9.]*' "$d/run_info.json" | awk '{print $2}')
  printf "%-18s %s%%\n" "$n" "$p"
done
