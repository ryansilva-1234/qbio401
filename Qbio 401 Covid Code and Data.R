# Load required libraries
library(DESeq2)        # differential expression analysis
library(pheatmap)      # heatmap visualization
library(EnhancedVolcano) # volcano plot
library(ggplot2)       # plotting framework
library(viridis)       # color scale for heatmap
library(RColorBrewer)  # palette for annotations

# Define relative data directory
data_dir <- file.path("data", "GSE171524_RAW")

# Import and combine count files
count_files <- list.files(
  path = data_dir,
  pattern = "_raw_counts\\.csv$",
  full.names = TRUE
)
cts_list <- lapply(count_files, function(f) {
  df <- read.csv(f, row.names = 1)
  setNames(df[, 2, drop = FALSE], basename(f))
})
counts <- do.call(cbind, cts_list)

# Build sample metadata
sample_names <- colnames(counts)
condition <- ifelse(grepl("ctr", sample_names), "control", "cov")
colData <- data.frame(
  row.names = sample_names,
  condition = factor(condition, levels = c("control", "cov"))
)

# Create DESeq2 dataset and filter low counts
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = colData,
  design    = ~ condition
)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "cov", "control"))
res_shrunk <- lfcShrink(dds, coef = "condition_cov_vs_control", type = "apeglm")
resOrdered <- res_shrunk[order(res_shrunk$padj), ]

# Save results
write.csv(as.data.frame(resOrdered), file = "DESeq2_results_cov_vs_ctrl.csv")

# Variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# PCA plot
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  stat_ellipse(aes(fill = condition), geom = "polygon", level = 0.95, alpha = 0.2, show.legend = FALSE) +
  geom_point(size = 4, stroke = 0.8) +
  scale_color_manual(values = c(control = "#66C2A5", cov = "#FC8D62")) +
  scale_fill_manual(values = c(control = "#66C2A5", cov = "#FC8D62")) +
  scale_shape_manual(values = c(control = 16, cov = 17)) +
  coord_fixed() +
  labs(
    title    = "PCA of Samples by Condition",
    subtitle = paste0("PC1: ", percentVar[1], "%   PC2: ", percentVar[2], "%"),
    x        = bquote(PC1~" ("~.(percentVar[1])~"%)"),
    y        = bquote(PC2~" ("~.(percentVar[2])~"%)")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle = element_text(size = 12, color = "grey30", hjust = 0),
    axis.title    = element_text(face = "bold", size = 14),
    axis.text     = element_text(size = 12)
  )

# Heatmap of top 20 genes by adjusted p-value
top20 <- head(order(res_shrunk$padj), 20)
mat_top20 <- assay(vsd)[top20, ]
orig_names  <- colnames(mat_top20)
short_names <- sub("GSM\\d+_([^_]+)_raw_counts\\.csv", "\\1", orig_names)
colnames(mat_top20) <- short_names

ann_col <- data.frame(
  Condition = colData[orig_names, "condition"]
)
rownames(ann_col) <- short_names
ord <- order(ann_col$Condition)
mat_top20 <- mat_top20[, ord]
ann_col   <- ann_col[ord, , drop = FALSE]

row_clusters <- cutree(hclust(dist(mat_top20)), k = 4)
annotation_row <- data.frame(Cluster = factor(row_clusters))
rownames(annotation_row) <- rownames(mat_top20)

ann_colors <- list(
  Condition = c(control = "skyblue", cov = "tomato"),
  Cluster   = setNames(brewer.pal(4, "Set2"), as.character(1:4))
)

pheatmap(
  mat_top20,
  scale             = "row",
  color             = viridis(100, option = "D"),
  annotation_col    = ann_col,
  annotation_row    = annotation_row,
  annotation_colors = ann_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  main              = "Top 20 DE Genes (z-score)"
)

# Volcano plot of differential expression
res_clean <- res[!is.na(res$pvalue) & !is.na(res$log2FoldChange), ]
top10_raw <- rownames(res_clean)[head(order(res_clean$pvalue), 10)]
EnhancedVolcano(
  res_shrunk,
  lab           = rownames(res_shrunk),
  x             = 'log2FoldChange',
  y             = 'pvalue',
  pCutoff       = 0.05,
  FCcutoff      = 1,
  selectLab     = top10_raw,
  title         = 'Volcano: cov vs control',
  subtitle      = 'p-value < 0.05 & |log₂FC| > 1',
  caption       = sprintf(
    'Total = %d; Up = %d; Down = %d',
    nrow(res_shrunk),
    sum(res_shrunk$pvalue < 0.05 & res_shrunk$log2FoldChange > 1, na.rm = TRUE),
    sum(res_shrunk$pvalue < 0.05 & res_shrunk$log2FoldChange < -1, na.rm = TRUE)
  )
) + theme_classic()