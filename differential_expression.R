#!/usr/bin/env Rscript
################################################################################
# RNA-seq Differential Expression Analysis using DESeq2
# 
# Author: Camila Avelino
# Date: 2024
# Description: Comprehensive differential expression analysis pipeline, including
#              normalization, quality control, statistical testing, and 
#              visualization of RNA-seq data
################################################################################

# Load required libraries --------------------------------------------------
library(tidyverse)      # Data manipulation and visualization
library(RColorBrewer)   # Color palettes
library(DESeq2)         # Differential expression analysis
library(pheatmap)       # Heatmaps
library(tximport)       # Import transcript-level estimates
library(apeglm)         # Log fold-change shrinkage
library(ggrepel)        # Better label placement in plots
library(GGally)         # Extended ggplot2 functionality

# Set up directories -------------------------------------------------------
base_dir <- "/path/to/your/project"  # MODIFY THIS
data_dir <- file.path(base_dir, "data/processed")
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(results_dir, "figures")
tables_dir <- file.path(results_dir, "tables")

# Create output directories if they don't exist
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# Parameters ---------------------------------------------------------------
ALPHA <- 0.05           # Significance threshold for p-values
S_VALUE_THRESHOLD <- 0.005  # Threshold for s-values (shrunken p-values)
LOG2FC_THRESHOLD <- 0   # Minimum log2 fold-change (0 = any change)

# ============================================================================
# PART 1: DATA IMPORT
# ============================================================================

message("Step 1: Importing data with tximport...")

# Load sample metadata
samples <- read.table(file.path(data_dir, "samples.txt"), 
                      header = TRUE, 
                      stringsAsFactors = FALSE)

# Construct file paths to RSEM output
files <- file.path(data_dir, 
                   paste0(samples$run, "_modificado.genes.results"))

# Name files for clarity
names(files) <- samples$sample_name

# Verify all files exist
if (!all(file.exists(files))) {
  missing <- files[!file.exists(files)]
  stop("Missing input files:\n", paste(missing, collapse = "\n"))
}

message("  Found ", length(files), " sample files")

# Import data using tximport
# This handles gene-level summarization and length scaling
txi <- tximport(files, 
                type = "rsem", 
                txIn = FALSE,  # Already at gene level
                txOut = FALSE)

# Load full metadata
meta <- read.table(file.path(data_dir, "meta.txt"), 
                   header = TRUE, 
                   row.names = 1,
                   stringsAsFactors = FALSE)

# Convert condition to factor with control as reference
meta$Condition <- factor(meta$Condition, 
                         levels = c("Control", "Patient"))

# Verify sample names match between count data and metadata
if (!all(colnames(txi$counts) == rownames(meta))) {
  stop("Sample names don't match between count matrix and metadata!")
}

message("  Successfully imported data for ", nrow(meta), " samples")
message("  Analyzing ", nrow(txi$counts), " genes")

# ============================================================================
# PART 2: CREATE DESEQ2 OBJECT AND NORMALIZATION
# ============================================================================

message("\nStep 2: Creating DESeq2 object and normalizing...")

# Create DESeqDataSet
dds <- DESeqDataSetFromTximport(txi = txi, 
                                colData = meta, 
                                design = ~ Condition)

# Pre-filtering: remove genes with very low counts
# Keep genes with at least 10 counts total across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

message("  After pre-filtering: ", nrow(dds), " genes retained")

# Estimate size factors for normalization
dds <- estimateSizeFactors(dds)

message("  Size factors calculated:")
print(sizeFactors(dds))

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Save normalized counts
write.csv(normalized_counts, 
          file = file.path(tables_dir, "normalized_counts.csv"))

# ============================================================================
# PART 3: QUALITY CONTROL
# ============================================================================

message("\nStep 3: Performing quality control...")

# Regularized log transformation for visualization
# blind=TRUE: transformation blind to experimental design
rld <- rlog(dds, blind = TRUE)
rld_mat <- assay(rld)

# --- PCA Plot ---
message("  Generating PCA plot...")

pca_data <- plotPCA(rld, 
                    intgroup = "Condition", 
                    ntop = nrow(rld),  # Use all genes
                    returnData = TRUE)

percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = name), size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("PCA - Sample Clustering")

ggsave(file.path(figures_dir, "PCA_plot.png"), 
       pca_plot, 
       width = 8, 
       height = 6,
       dpi = 300)

# --- Sample Correlation Heatmap ---
message("  Generating correlation heatmap...")

rld_cor <- cor(rld_mat)

png(file.path(figures_dir, "correlation_heatmap.png"), 
    width = 10, 
    height = 10, 
    units = "in", 
    res = 300)

pheatmap(rld_cor,
         display_numbers = TRUE,
         number_color = "black",
         fontsize_number = 8,
         main = "Sample-to-Sample Correlation",
         border_color = NA)

dev.off()

# ============================================================================
# PART 4: DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================

message("\nStep 4: Running differential expression analysis...")

# Run DESeq2 analysis
# This performs: estimation of size factors, dispersion, and statistical testing
dds <- DESeq(dds)

# Check results names
message("  Available contrasts:")
print(resultsNames(dds))

# Plot dispersion estimates
png(file.path(figures_dir, "dispersion_plot.png"), 
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)

plotDispEsts(dds, 
             main = "Dispersion Estimates")

dev.off()

# Extract results
# Note: DESeq2 automatically uses "Control" as reference level
results_unshrunken <- results(dds, 
                              contrast = c("Condition", "Patient", "Control"),
                              alpha = ALPHA)

# Apply log fold-change shrinkage (apeglm)
# This reduces noise in lowly expressed genes
results_shrunken <- lfcShrink(dds, 
                              coef = "Condition_Patient_vs_Control",
                              type = "apeglm",
                              svalue = TRUE)  # Calculate s-values

# Order by adjusted p-value
results_shrunken <- results_shrunken[order(results_shrunken$svalue), ]

# Summary
message("\n  === RESULTS SUMMARY ===")
summary(results_shrunken, alpha = S_VALUE_THRESHOLD)

# Save complete results
write.table(as.data.frame(results_shrunken),
            file = file.path(tables_dir, "deseq2_results_all.txt"),
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# --- MA Plot ---
message("  Generating MA plot...")

png(file.path(figures_dir, "MA_plot.png"), 
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)

plotMA(results_shrunken, 
       ylim = c(-5, 5),
       alpha = S_VALUE_THRESHOLD,
       main = "MA Plot: Patient vs Control")

dev.off()

# ============================================================================
# PART 5: EXTRACT SIGNIFICANT GENES
# ============================================================================

message("\nStep 5: Extracting significant genes...")

# Remove genes with NA values
results_complete <- na.omit(results_shrunken)

message("  Genes after removing NAs: ", nrow(results_complete))

# Define significant genes
sig_genes <- results_complete[results_complete$svalue < S_VALUE_THRESHOLD, ]

message("  Total significant genes: ", nrow(sig_genes))

# Separate upregulated and downregulated
upregulated <- sig_genes[sig_genes$log2FoldChange > LOG2FC_THRESHOLD, ]
downregulated <- sig_genes[sig_genes$log2FoldChange < -LOG2FC_THRESHOLD, ]

message("  Upregulated in patients: ", nrow(upregulated))
message("  Downregulated in patients: ", nrow(downregulated))

# Save significant gene lists
write.table(as.data.frame(upregulated),
            file = file.path(tables_dir, "upregulated_genes.txt"),
            sep = "\t",
            quote = FALSE,
            col.names = NA)

write.table(as.data.frame(downregulated),
            file = file.path(tables_dir, "downregulated_genes.txt"),
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Create a combined list with direction labels
sig_genes_annotated <- as.data.frame(sig_genes)
sig_genes_annotated$direction <- ifelse(sig_genes_annotated$log2FoldChange > 0,
                                        "Upregulated",
                                        "Downregulated")

write.csv(sig_genes_annotated,
          file = file.path(tables_dir, "significant_genes_annotated.csv"))

# ============================================================================
# PART 6: VISUALIZATION - HEATMAP
# ============================================================================

message("\nStep 6: Generating expression heatmap...")

# Extract normalized expression for significant genes
norm_sig <- normalized_counts[rownames(sig_genes), ]

# Create annotation for heatmap
annotation_col <- data.frame(
  Condition = meta$Condition,
  row.names = rownames(meta)
)

# Set color palette
heat_colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

# Generate heatmap
png(file.path(figures_dir, "expression_heatmap.png"), 
    width = 10, 
    height = 12, 
    units = "in", 
    res = 300)

pheatmap(norm_sig,
         color = heat_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,  # Too many genes to show
         annotation_col = annotation_col,
         border_color = NA,
         fontsize = 10,
         scale = "row",  # Z-score scaling by row
         main = paste0("Expression Heatmap\n(", 
                       nrow(sig_genes), 
                       " significant genes)"))

dev.off()

# --- Top 50 genes heatmap ---
top50 <- head(sig_genes, 50)
norm_top50 <- normalized_counts[rownames(top50), ]

png(file.path(figures_dir, "expression_heatmap_top50.png"), 
    width = 10, 
    height = 10, 
    units = "in", 
    res = 300)

pheatmap(norm_top50,
         color = heat_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = annotation_col,
         border_color = NA,
         fontsize = 10,
         fontsize_row = 8,
         scale = "row",
         main = "Top 50 Most Significant Genes")

dev.off()

# ============================================================================
# PART 7: VOLCANO PLOT
# ============================================================================

message("\nStep 7: Generating volcano plot...")

# Prepare data for volcano plot
volcano_data <- as.data.frame(results_complete)
volcano_data$significant <- volcano_data$svalue < S_VALUE_THRESHOLD
volcano_data$gene <- rownames(volcano_data)

# Label top genes
volcano_data$label <- ""
top_genes <- head(volcano_data[order(volcano_data$svalue), ], 20)
volcano_data$label[rownames(volcano_data) %in% rownames(top_genes)] <- 
  rownames(top_genes)

# Create volcano plot
volcano_plot <- ggplot(volcano_data, 
                       aes(x = log2FoldChange, 
                           y = -log10(svalue),
                           color = significant)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(aes(label = label), 
                  size = 3,
                  max.overlaps = 20) +
  geom_hline(yintercept = -log10(S_VALUE_THRESHOLD), 
             linetype = "dashed",
             color = "blue") +
  geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD),
             linetype = "dashed",
             color = "blue") +
  labs(title = "Volcano Plot: Patient vs Control",
       x = "Log2 Fold Change",
       y = "-Log10 S-value") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "volcano_plot.png"),
       volcano_plot,
       width = 10,
       height = 8,
       dpi = 300)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

message("\n============================================")
message("ANALYSIS COMPLETE!")
message("============================================")
message("\nResults saved to: ", results_dir)
message("\nSummary Statistics:")
message("  Total genes analyzed: ", nrow(results_complete))
message("  Significant genes (s-value < ", S_VALUE_THRESHOLD, "): ", nrow(sig_genes))
message("    - Upregulated: ", nrow(upregulated))
message("    - Downregulated: ", nrow(downregulated))
message("\nOutput files:")
message("  - Figures: ", figures_dir)
message("  - Tables: ", tables_dir)
message("\n============================================")

# Save session info for reproducibility
sink(file.path(results_dir, "session_info.txt"))
sessionInfo()
sink()

message("\nSession info saved to: ", file.path(results_dir, "session_info.txt"))
