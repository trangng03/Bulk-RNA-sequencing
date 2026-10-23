# Paired-end RNA-seq pipeline

This template runs **FastQC → STAR alignment and gene counts → pause to check strandedness → DESeq2 → fgsea**. It expects gzipped paired-end reads, a prebuilt STAR genome index made with the matching GTF, and a gene-symbol GMT gene-set file. It does not download references or install software.

## Before running

1. Activate your existing `rnaseq` conda environment: `conda activate rnaseq`.
2. Edit `samples.tsv` to give **real absolute FASTQ paths** and correct biological conditions. The first three sample names follow your control1–control3 naming; `treated1`–`treated3` and their labels are placeholders until you confirm the experimental groups. Use a tab-separated file with four columns exactly.
3. Provide the path to your existing STAR genome index, its matching human GTF, and a human gene-symbol GMT file (for example an MSigDB Hallmark symbols GMT acquired under its terms). Keep the genome build and annotation consistent. The same gene symbols must appear in the GTF and GMT.
4. After STAR, confirm library strandedness. Use `0` for unstranded, `1` for forward, `2` for reverse. STAR's `ReadsPerGene.out.tab` columns 2/3/4 correspond respectively. The script prints totals for each sample, then pauses. Compare at least one sample from each condition; when forward/reverse totals are ambiguous, confirm with library-preparation records or a strand-inference tool before DESeq2.
5. Check the sample sheet reflects the study design (e.g. batch or pairing). This script fits `~ Condition` only. Adjust that design before interpreting results if other factors need accounting for.

## Run

From this directory, first align the reads (replace the index path):

```bash
conda activate rnaseq
bash run_pipeline.sh samples.tsv /absolute/path/star_index 8
```

Inspect the totals printed at the end and choose the appropriate setting. For instance, if the reverse column dominates consistently (`2`):

```bash
bash run_analysis.sh results /absolute/path/annotation.gtf /absolute/path/h.all.symbols.gmt treated control 2
```

The `2` is only an example. If forward/reverse totals are similar and the library is unstranded, use `0`; if forward dominates, use `1`. By default results go to `./results`. To choose another directory: `OUT_DIR=/absolute/path/project/results bash run_pipeline.sh ...`, then pass that same directory as the first argument to `run_analysis.sh`.

The scripts stop and report missing commands (`fastqc`, `STAR`, `Rscript`) or missing R packages (`DESeq2`, `fgsea`). Share that error before installing anything. Inspect `results/logs/` for per-sample logs and `results/02_star/<sample>/Log.final.out` for mapping quality.

## Folder layout

| Folder | Files |
| --- | --- |
| `01_fastqc/<sample>/` | Per-read FastQC HTML and ZIP reports |
| `02_star/<sample>/` | Sorted BAM, STAR gene counts, alignment logs |
| `03_counts/` | Input manifest and combined raw counts |
| `04_deseq2/` | Full/significant differential expression tables, normalized counts, PCA and MA plots, saved DESeq2 object |
| `05_gsea/` | Full ranked gene list and pathway enrichment results |
| `logs/` | FastQC, STAR, R analysis console output |

DESeq2 uses raw integer counts, not normalized counts. The reported log2 fold change is **treated over control** in the example command. Significance table uses adjusted p < 0.05 and |log2 fold change| ≥ 1; all tested genes are retained separately. GSEA ranks all genes with finite Wald statistics and recognized symbols, irrespective of differential-expression significance. The `leadingEdge` column lists core genes separated by semicolons. STAR GTF gene IDs are mapped to gene symbols by matching gene IDs in the annotation GTF; version suffixes are removed for this lookup.

This is a two-condition, unpaired design with at least two biological replicates per group. The workflow does not merge technical runs or run MultiQC. Rerunning `run_analysis.sh` with a corrected strand choice uses the existing STAR counts without repeating alignment. Rerunning `run_pipeline.sh` into the same output location repeats alignment and overwrites many earlier outputs; choose a new `OUT_DIR` for a new analysis.
