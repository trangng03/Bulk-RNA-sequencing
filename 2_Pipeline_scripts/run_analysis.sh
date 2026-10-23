#!/usr/bin/env bash
set -euo pipefail
# Usage: bash run_analysis.sh RESULTS_DIR annotation.gtf pathways.gmt CASE CONTROL STRAND
if (( $# != 6 )); then
  echo "Usage: bash $0 RESULTS_DIR annotation.gtf pathways.gmt CASE CONTROL STRAND" >&2
  exit 2
fi
OUT=$(realpath "$1"); GTF=$(realpath "$2"); GMT=$(realpath "$3")
CASE=$4; CONTROL=$5; STRAND=$6
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
command -v Rscript >/dev/null || { echo "Rscript missing in active environment; please tell me before installing." >&2; exit 1; }
Rscript -e 'x <- c("DESeq2", "fgsea"); missing <- x[!vapply(x, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop("Missing R packages: ", paste(missing, collapse=", "), ". Please ask before installing anything.")'
for f in "$OUT/03_counts/input_manifest.tsv" "$GTF" "$GMT"; do [[ -f "$f" ]] || { echo "File missing: $f" >&2; exit 1; }; done
[[ "$STRAND" =~ ^[012]$ ]] || { echo "STRAND must be 0, 1, or 2" >&2; exit 1; }
[[ "$CASE" != "$CONTROL" ]] || { echo "CASE and CONTROL must differ" >&2; exit 1; }
Rscript "$SCRIPT_DIR/analyze.R" "$OUT/03_counts/input_manifest.tsv" "$OUT" "$GTF" "$GMT" "$CASE" "$CONTROL" "$STRAND" \
  > "$OUT/logs/analysis.log" 2>&1 || { cat "$OUT/logs/analysis.log" >&2; exit 1; }
echo "DESeq2 and GSEA done: $OUT"
