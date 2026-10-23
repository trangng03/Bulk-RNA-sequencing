#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop("Expected: manifest.tsv results_dir annotation.gtf pathways.gmt CASE CONTROL STRAND")
manifest <- args[1]; out <- args[2]; gtf <- args[3]; gmt <- args[4]
case <- args[5]; control <- args[6]; strand <- as.integer(args[7])
if (!strand %in% 0:2) stop("STRAND must be 0, 1, or 2")
samples <- read.delim(manifest, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("Sample", "Condition") %in% names(samples)), !anyDuplicated(samples$Sample)) 
if (!all(c(case, control) %in% samples$Condition)) stop("CASE/CONTROL absent from sample sheet")
if (any(!samples$Condition %in% c(case, control))) stop("This script analyzes exactly two conditions; check samples.tsv")
if (any(table(samples$Condition) < 2)) stop("Use at least two biological replicates in each condition")

paths <- file.path(out, "02_star", samples$Sample, "ReadsPerGene.out.tab")
if (!all(file.exists(paths))) stop("Missing STAR gene counts: ", paste(paths[!file.exists(paths)], collapse=", "))
counts <- lapply(paths, function(path) {

  x <- read.delim(path, header=FALSE, check.names=FALSE, stringsAsFactors=FALSE)
  if (ncol(x) != 4L || nrow(x) < 5L) stop("Unexpected STAR GeneCounts file: ", path)
  x <- x[-(1:4), , drop=FALSE]  # Four STAR summary rows are not genes.
  v <- as.numeric(x[[strand + 2L]])
  if (anyNA(v) || any(v < 0)) stop("Invalid counts in ", path)
  setNames(v, x[[1L]])
})
ids <- names(counts[[1L]])
if (any(vapply(counts, function(x) !identical(names(x), ids), logical(1)))) stop("Gene IDs/order differ across STAR files; check that the same index was used")
cts <- do.call(cbind, counts)
rownames(cts) <- ids; colnames(cts) <- samples$Sample
if (anyDuplicated(ids)) stop("Duplicate gene IDs in STAR output")
write.csv(data.frame(gene_id=ids, cts, check.names=FALSE), file.path(out,"03_counts","raw_gene_counts.csv"), row.names=FALSE)

samples$Condition <- factor(samples$Condition, levels=c(control, case))
rownames(samples) <- samples$Sample
dds <- DESeq2::DESeqDataSetFromMatrix(countData=cts, colData=samples, design=~Condition)
dds <- dds[rowSums(DESeq2::counts(dds) >= 10) >= 2, ]
if (nrow(dds) == 0L) stop("No genes passed count filter; check strandedness and STAR logs")
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast=c("Condition", case, control))
resdf <- data.frame(gene_id=rownames(res), as.data.frame(res), check.names=FALSE)
resdf <- resdf[order(resdf$padj, na.last=TRUE), , drop=FALSE]

write.csv(data.frame(gene_id=rownames(dds), DESeq2::counts(dds, normalized=TRUE), check.names=FALSE),
          file.path(out,"04_deseq2","normalized_counts.csv"), row.names=FALSE)
saveRDS(dds, file.path(out,"04_deseq2","deseq2_object.rds"))
pdf(file.path(out,"04_deseq2","MA_plot.pdf")); DESeq2::plotMA(res, ylim=c(-5,5)); dev.off()
vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)
pca <- prcomp(t(SummarizedExperiment::assay(vsd)))
write.csv(data.frame(sample=samples$Sample, condition=samples$Condition, PC1=pca$x[,1], PC2=pca$x[,2]),
          file.path(out,"04_deseq2","PCA_coordinates.csv"), row.names=FALSE)
pdf(file.path(out,"04_deseq2","PCA.pdf"))
plot(pca$x[,1:2], col=as.integer(samples$Condition), pch=19, xlab="PC1", ylab="PC2")
text(pca$x[,1:2], labels=samples$Sample, pos=3, cex=.7)
legend("topright", legend=levels(samples$Condition), col=seq_along(levels(samples$Condition)), pch=19)
dev.off()

# Read only gene features from the GTF, in chunks, to map Ensembl IDs to gene symbols.
con <- if (grepl("\\.gz$", gtf)) gzfile(gtf, "rt") else file(gtf, "rt")
gene_id <- character(); gene_symbol <- character()
extract <- function(z, key) {
  m <- regexec(paste0('(?:^|;[[:space:]]*)', key, '[[:space:]]+"([^"]+)"'), z, perl=TRUE)
  vapply(regmatches(z, m), function(a) if (length(a) > 1L) a[2L] else NA_character_, character(1))
}
repeat {
  lines <- readLines(con, n=50000L, warn=FALSE)
  if (!length(lines)) break
  lines <- lines[grepl("^[^#].*\\tgene\\t", lines)]
  if (!length(lines)) next
  attrs <- sub("^(?:[^\\t]*\\t){8}", "", lines, perl=TRUE)
  gene_id <- c(gene_id, extract(attrs, "gene_id"))
  gene_symbol <- c(gene_symbol, extract(attrs, "gene_name"))
}
close(con)
gene_map <- setNames(gene_symbol, sub("\\.[0-9]+$", "", gene_id))
gene_map <- gene_map[!duplicated(names(gene_map))]
gene_names <- unname(gene_map[sub("\\.[0-9]+$", "", resdf$gene_id)])
missing_names <- is.na(gene_names) | !nzchar(gene_names)
gene_names[missing_names] <- resdf$gene_id[missing_names]
output_df <- resdf
output_df$gene_id <- gene_names
names(output_df)[names(output_df) == "gene_id"] <- "gene_name"
write.csv(output_df, file.path(out, "04_deseq2", "all_genes.csv"), row.names = FALSE)
write.csv(subset(output_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >=1),
          file.path(out, "04_deseq2", "significant_genes.csv"), row.names = FALSE)
# ---- Volcano plot ----
needed <- c("ggplot2", "ggrepel")
missing_packages <- needed[
  !vapply(needed, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

# Use the results containing gene names.
volcano_df <- output_df

# Keep rows with usable fold changes and adjusted p-values.
volcano_df <- volcano_df[
  is.finite(volcano_df$log2FoldChange) &
    is.finite(volcano_df$padj) &
    volcano_df$padj >= 0,
  ,
  drop = FALSE
]

if (nrow(volcano_df) == 0) {
  stop("No genes have usable values for the volcano plot")
}

# Calculate the y-axis; avoid infinity when padj is zero.
volcano_df$neg_log10_padj <- -log10(
  pmax(volcano_df$padj, .Machine$double.xmin)
)

# Classify genes using your existing thresholds.
volcano_df$direction <- "Below thresholds"

volcano_df$direction[
  volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange >= 1
] <- "Upregulated"

volcano_df$direction[
  volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange <= -1
] <- "Downregulated"

# Find the top 5 significant genes in each direction.
up <- subset(volcano_df, direction == "Upregulated")
down <- subset(volcano_df, direction == "Downregulated")

up <- up[order(-up$log2FoldChange, up$padj), , drop = FALSE]
down <- down[order(down$log2FoldChange, down$padj), , drop = FALSE]

top_up <- head(up, 5)
top_down <- head(down, 5)
label_df <- rbind(top_up, top_down)

# Assign labels only to the selected rows.
volcano_df$label <- ""
selected <- rownames(volcano_df) %in% rownames(label_df)
volcano_df$label[selected] <- volcano_df$gene_name[selected]

# Create the plot.
volcano_plot <- ggplot2::ggplot(
  volcano_df,
  ggplot2::aes(x = log2FoldChange, y = neg_log10_padj)
) +
  ggplot2::geom_point(
    ggplot2::aes(color = direction),
    alpha = 0.65,
    size = 1.6
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Upregulated" = "#D73027",
      "Downregulated" = "#4575B4",
      "Below thresholds" = "grey75"
    )
  ) +
  ggplot2::geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey50"
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50"
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = label),
    seed = 42,
    max.overlaps = Inf,
    box.padding = 0.5,
    min.segment.length = 0,
    size = 3.5
  ) +
  ggplot2::labs(
    title = paste(case, "vs", control),
    subtitle = "Top 5 significant genes per direction by log2 fold change",
    x = "Log2 fold change",
    y = expression(-log[10]("adjusted p-value")),
    color = NULL
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "bottom")

# Save PNG and PDF.
ggplot2::ggsave(
  file.path(out, "04_deseq2", "volcano_plot.png"),
  plot = volcano_plot,
  width = 10, height = 8, dpi = 300
)

ggplot2::ggsave(
  file.path(out, "04_deseq2", "volcano_plot.pdf"),
  plot = volcano_plot,
  width = 10, height = 8
)

# Save the labeled genes.
write.csv(
  label_df[, c(
    "gene_name", "log2FoldChange", "padj", "direction"
  ), drop = FALSE],
  file.path(out, "04_deseq2", "volcano_top_genes.csv"),
  row.names = FALSE
)
rankdf <- resdf[is.finite(resdf$stat), c("gene_id", "stat")]
rankdf$symbol <- unname(gene_map[sub("\\.[0-9]+$", "", rankdf$gene_id)])
rankdf <- rankdf[!is.na(rankdf$symbol) & nzchar(rankdf$symbol), ]
rankdf <- rankdf[order(abs(rankdf$stat), decreasing=TRUE), ]
rankdf <- rankdf[!duplicated(rankdf$symbol), ]
if (nrow(rankdf) < 15L) stop("Too few gene symbols matched the GTF; check reference compatibility")
ranks <- setNames(rankdf$stat, rankdf$symbol)
ranks <- sort(ranks, decreasing=TRUE)
write.csv(data.frame(symbol=names(ranks), wald_stat=unname(ranks)),
          file.path(out,"05_gsea","ranked_genes.csv"), row.names=FALSE)
pathways <- fgsea::gmtPathways(gmt)
set.seed(42)
enrichment <- fgsea::fgsea(pathways=pathways, stats=ranks, minSize=15, maxSize=500)
enrichment <- as.data.frame(enrichment)
if (nrow(enrichment)) {
  enrichment$leadingEdge <- vapply(enrichment$leadingEdge, paste, collapse=";", FUN.VALUE=character(1))
  enrichment <- enrichment[order(enrichment$padj), ]
}
write.csv(enrichment, file.path(out,"05_gsea","enrichment.csv"), row.names=FALSE)
cat("Done:", nrow(dds), "tested genes;", length(ranks), "genes ranked by symbol;", nrow(enrichment), "pathways tested.\n")
