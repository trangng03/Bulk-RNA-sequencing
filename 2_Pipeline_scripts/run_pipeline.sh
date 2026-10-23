#!/usr/bin/env bash
set -euo pipefail

# Run inside your existing rnaseq environment. No packages are installed here.
# Usage: bash run_pipeline.sh samples.tsv STAR_INDEX [THREADS]
# Pause after STAR so strandedness can be checked before DESeq2.
if (( $# < 2 || $# > 3 )); then
  echo "Usage: bash $0 samples.tsv STAR_INDEX [THREADS]" >&2
  exit 2
fi
SAMPLES=$(realpath "$1")
INDEX=$(realpath "$2")
THREADS=${3:-8}
OUT=$(realpath -m "${OUT_DIR:-$PWD/results}")
for cmd in fastqc STAR realpath; do
  command -v "$cmd" >/dev/null || { echo "Missing $cmd in the active environment. Please tell me before installing anything." >&2; exit 1; }
done
[[ -f "$SAMPLES" ]] || { echo "File missing: $SAMPLES" >&2; exit 1; }
[[ -f "$INDEX/genomeParameters.txt" ]] || { echo "Not a STAR genome index: $INDEX" >&2; exit 1; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 1; }

mkdir -p "$OUT"/{01_fastqc,02_star,03_counts,04_deseq2,05_gsea,logs}
printf 'Sample\tCondition\tRead1\tRead2\n' > "$OUT/03_counts/input_manifest.tsv"
declare -A seen=()
num=0
while IFS=$'\t' read -r sample condition r1 r2 extra || [[ -n "${sample:-}" ]]; do
  [[ "$sample" == sample && "$condition" == condition ]] && continue
  [[ -z "${extra:-}" && "$sample" =~ ^[A-Za-z0-9_.-]+$ && -n "$condition" && -n "$r1" && -n "$r2" ]] || { echo "Invalid row in $SAMPLES: expected sample, condition, fastq_1, fastq_2 (tab separated); simple sample IDs only" >&2; exit 1; }
  [[ -z "${seen[$sample]:-}" ]] || { echo "Duplicate sample: $sample" >&2; exit 1; }
  seen[$sample]=1
  for fq in "$r1" "$r2"; do
    [[ -f "$fq" && "$fq" == *.gz ]] || { echo "Missing gzipped FASTQ: $fq (use absolute paths in samples.tsv)" >&2; exit 1; }
  done
  printf '%s\t%s\t%s\t%s\n' "$sample" "$condition" "$r1" "$r2" >> "$OUT/03_counts/input_manifest.tsv"
  ((num+=1))
done < "$SAMPLES"
(( num >= 2 )) || { echo "Need at least two samples" >&2; exit 1; }

echo "FastQC -> $OUT/01_fastqc"
tail -n +2 "$OUT/03_counts/input_manifest.tsv" | while IFS=$'\t' read -r sample condition r1 r2; do
  mkdir -p "$OUT/01_fastqc/$sample" "$OUT/02_star/$sample"
  fastqc -t 2 -o "$OUT/01_fastqc/$sample" "$r1" "$r2" > "$OUT/logs/${sample}.fastqc.log" 2>&1
  echo "STAR: $sample"
  STAR --runThreadN "$THREADS" --genomeDir "$INDEX" \
    --readFilesIn "$r1" "$r2" --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts \
    --outFileNamePrefix "$OUT/02_star/$sample/" \
    > "$OUT/logs/${sample}.star.log" 2>&1
done
echo "Alignment done. Compare column totals for at least one control and one treated sample:"
for sample in "${!seen[@]}"; do
  awk -v sample="$sample" 'NR>4 {a+=$2; b+=$3; c+=$4} END {printf "%s: unstranded=%.0f forward=%.0f reverse=%.0f\n", sample,a,b,c}' "$OUT/02_star/$sample/ReadsPerGene.out.tab"
done
echo "Choose STRAND=0 (unstranded), 1 (forward), or 2 (reverse), then run run_analysis.sh. Results: $OUT"
