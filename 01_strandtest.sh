#!/bin/bash
# Determines strandedness by quantifying one library three ways.
# Result: unstranded 92.3%, fr 1.0%, rf 88.8% -> reverse-stranded.
set -u
K="$HOME/Downloads/kallisto/kallisto"; IDX="$HOME/Downloads/index.idx"
OUT="$HOME/postn_rnaseq"; mkdir -p "$OUT"
R1=$(find "$HOME/Library/CloudStorage" -name "3818_Sample1_*_R1_001.fastq.gz" | head -1)
R2=$(find "$HOME/Library/CloudStorage" -name "3818_Sample1_*_R2_001.fastq.gz" | head -1)
for mode in un fr rf; do
  case "$mode" in un) F="" ;; fr) F="--fr-stranded" ;; rf) F="--rf-stranded" ;; esac
  "$K" quant -i "$IDX" -o "$OUT/st_$mode" -t 4 $F "$R1" "$R2"
done
for mode in un fr rf; do
  printf "%-3s " "$mode"; grep p_pseudoaligned "$OUT/st_$mode/run_info.json"
done
