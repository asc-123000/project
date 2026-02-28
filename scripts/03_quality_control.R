################################################################################
# 03_quality_control.R - Comprehensive Quality Control Analysis
################################################################################
#
# Project: TP4-Induced Transcriptomic Changes in Triple-Negative Breast Cancer
# Dataset: GSE74764 (Agilent SurePrint G3 Human GE v2 8x60K)
# 
# Purpose:
#   Perform thorough quality control to:
#   1. Verify data quality and normalization
#   2. Assess sample relationships (PCA, clustering)
#   3. Identify potential outliers
#   4. Confirm experimental design is reflected in data structure
#   5. Generate publication-quality QC figures
#
# Scientific Rationale:
# ─────────────────────────────────────────────────
#   QC is critical for microarray analysis because:
#   
#   1. TECHNICAL VARIATION: Batch effects, hybridization issues, or
#      degraded samples can obscure biological signals.
#   
#   2. SAMPLE VERIFICATION: PCA should show samples grouping by
#      biological factors (cell line, treatment), not technical artifacts.
#   
#   3. SMALL SAMPLE SIZE: With n=2 per group, a single outlier can
#      dramatically affect results. Early detection is essential.
#   
#   4. NORMALIZATION VERIFICATION: We're using pre-processed data,
#      so we must verify the normalization is adequate.
#
# Author: Computational Biology Pipeline
# Date: 2026-01-19
#
################################################################################

# =============================================================================
# 1. SETUP
# =============================================================================

# Source the setup script
source("scripts/00_setup.R")

# Load required packages
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(cowplot)
library(tidyverse)

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Script 03: Quality Control Analysis                                     ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

# =============================================================================
# 2. LOAD PREPROCESSED DATA
# =============================================================================

message("--- Loading Preprocessed Data ---\n")

# Load preprocessing output
preproc_output <- readRDS(file.path(PATHS$processed, "02_preprocessing_output.rds"))

expr_probe <- preproc_output$expression_probe
expr_gene <- preproc_output$expression_gene
sample_metadata <- preproc_output$sample_metadata
sample_stats <- preproc_output$sample_stats

# Use gene-level expression if available, otherwise probe-level
if (!is.null(expr_gene)) {
  expr_matrix <- expr_gene
  message("Using gene-level expression matrix: ", nrow(expr_matrix), " genes")
} else {
  expr_matrix <- expr_probe
  message("Using probe-level expression matrix: ", nrow(expr_matrix), " probes")
}

message("Samples: ", ncol(expr_matrix))
message("Groups: ", paste(levels(sample_metadata$group), collapse = ", "))

# =============================================================================
# 3. EXPRESSION DISTRIBUTION ANALYSIS
# =============================================================================

message("\n--- Expression Distribution Analysis ---\n")

# 3.1 Density plots
message("Creating density plots...")

# Prepare data for ggplot
expr_long <- expr_matrix %>%
  as.data.frame() %>%
  mutate(feature = rownames(.)) %>%
  pivot_longer(-feature, names_to = "sample_id", values_to = "expression") %>%
  left_join(sample_metadata[, c("sample_id", "cell_line", "treatment", "group")], 
            by = "sample_id")

# Density plot colored by sample
p_density_sample <- ggplot(expr_long, aes(x = expression, color = sample_id)) +
  geom_density(linewidth = 0.8, alpha = 0.8) +
  labs(
    title = "Expression Density Distribution by Sample",
    x = "Log2 Expression",
    y = "Density",
    color = "Sample"
  ) +
  theme_publication() +
  theme(legend.position = "right")

# Density plot colored by cell line
p_density_cellline <- ggplot(expr_long, aes(x = expression, color = cell_line)) +
  geom_density(linewidth = 1, alpha = 0.8) +
  scale_color_manual(values = ANALYSIS_PARAMS$colors$cell_line) +
  labs(
    title = "Expression Density Distribution by Cell Line",
    x = "Log2 Expression",
    y = "Density",
    color = "Cell Line"
  ) +
  theme_publication()

# Density plot colored by treatment
p_density_treatment <- ggplot(expr_long, aes(x = expression, color = treatment)) +
  geom_density(linewidth = 1, alpha = 0.8) +
  scale_color_manual(values = ANALYSIS_PARAMS$colors$treatment) +
  labs(
    title = "Expression Density Distribution by Treatment",
    x = "Log2 Expression",
    y = "Density",
    color = "Treatment"
  ) +
  theme_publication()

# Density plot by group (faceted)
p_density_group <- ggplot(expr_long, aes(x = expression, fill = group)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~group, ncol = 2) +
  labs(
    title = "Expression Density Distribution by Group",
    x = "Log2 Expression",
    y = "Density"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# Save density plots
ggsave(file.path(PATHS$figures_qc, "density_by_sample.pdf"), 
       p_density_sample, width = 10, height = 6, dpi = ANALYSIS_PARAMS$plot_dpi)
ggsave(file.path(PATHS$figures_qc, "density_by_sample.png"), 
       p_density_sample, width = 10, height = 6, dpi = ANALYSIS_PARAMS$plot_dpi)

ggsave(file.path(PATHS$figures_qc, "density_by_cellline.pdf"), 
       p_density_cellline, width = 8, height = 6, dpi = ANALYSIS_PARAMS$plot_dpi)

ggsave(file.path(PATHS$figures_qc, "density_by_treatment.pdf"), 
       p_density_treatment, width = 8, height = 6, dpi = ANALYSIS_PARAMS$plot_dpi)

message("✓ Density plots saved")

# 3.2 Boxplots
message("Creating boxplots...")

# Boxplot of expression by sample
p_box_sample <- ggplot(expr_long, aes(x = sample_id, y = expression, fill = group)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  labs(
    title = "Expression Distribution by Sample",
    x = "Sample",
    y = "Log2 Expression",
    fill = "Group"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Boxplot by group
p_box_group <- ggplot(expr_long, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  labs(
    title = "Expression Distribution by Group",
    x = "Group",
    y = "Log2 Expression"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Save boxplots
ggsave(file.path(PATHS$figures_qc, "boxplot_by_sample.pdf"), 
       p_box_sample, width = 10, height = 6, dpi = ANALYSIS_PARAMS$plot_dpi)
ggsave(file.path(PATHS$figures_qc, "boxplot_by_sample.png"), 
       p_box_sample, width = 10, height = 6, dpi = ANALYSIS_PARAMS$plot_dpi)

ggsave(file.path(PATHS$figures_qc, "boxplot_by_group.pdf"), 
       p_box_group, width = 8, height = 6, dpi = ANALYSIS_PARAMS$plot_dpi)

message("✓ Boxplots saved")

# =============================================================================
# 4. PRINCIPAL COMPONENT ANALYSIS
# =============================================================================

message("\n--- Principal Component Analysis ---\n")

# Perform PCA on transposed expression matrix (samples as rows)
pca_result <- prcomp(t(expr_matrix), center = TRUE, scale. = TRUE)

# Calculate variance explained
var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

message("Variance explained by first 5 PCs:")
for (i in 1:min(5, length(var_explained))) {
  message("  PC", i, ": ", round(var_explained[i], 1), "%")
}

# Create PCA data frame for plotting
pca_df <- data.frame(
  sample_id = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  PC4 = pca_result$x[, 4]
)

# Add metadata
pca_df <- merge(pca_df, sample_metadata, by = "sample_id")

# 4.1 PCA plot - PC1 vs PC2 (colored by cell line, shaped by treatment)
p_pca_main <- ggplot(pca_df, aes(x = PC1, y = PC2, 
                                   color = cell_line, shape = treatment)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 20) +
  scale_color_manual(values = ANALYSIS_PARAMS$colors$cell_line) +
  scale_shape_manual(values = c("Mock" = 16, "TP4" = 17)) +
  labs(
    title = "PCA of Gene Expression",
    subtitle = paste0("GSE74764: TP4 treatment in TNBC vs Fibroblasts"),
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)"),
    color = "Cell Line",
    shape = "Treatment"
  ) +
  theme_publication() +
  theme(
    legend.position = "right",
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

# 4.2 PCA plot - colored by group (single factor)
p_pca_group <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 20) +
  labs(
    title = "PCA of Gene Expression by Group",
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)"),
    color = "Group"
  ) +
  theme_publication()

# 4.3 PCA plot - PC3 vs PC4
p_pca_34 <- ggplot(pca_df, aes(x = PC3, y = PC4, 
                                 color = cell_line, shape = treatment)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 20) +
  scale_color_manual(values = ANALYSIS_PARAMS$colors$cell_line) +
  scale_shape_manual(values = c("Mock" = 16, "TP4" = 17)) +
  labs(
    title = "PCA - Higher Components",
    x = paste0("PC3 (", round(var_explained[3], 1), "% variance)"),
    y = paste0("PC4 (", round(var_explained[4], 1), "% variance)"),
    color = "Cell Line",
    shape = "Treatment"
  ) +
  theme_publication()

# 4.4 Scree plot
scree_df <- data.frame(
  PC = factor(paste0("PC", 1:length(var_explained)), 
              levels = paste0("PC", 1:length(var_explained))),
  Variance = var_explained,
  Cumulative = cumsum(var_explained)
)

p_scree <- ggplot(scree_df[1:8, ], aes(x = PC)) +
  geom_bar(aes(y = Variance), stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_line(aes(y = Cumulative, group = 1), color = "red", linewidth = 1) +
  geom_point(aes(y = Cumulative), color = "red", size = 3) +
  scale_y_continuous(
    name = "Variance Explained (%)",
    sec.axis = sec_axis(~., name = "Cumulative Variance (%)")
  ) +
  labs(
    title = "PCA Scree Plot",
    x = "Principal Component"
  ) +
  theme_publication()

# Save PCA plots
ggsave(file.path(PATHS$figures_qc, "PCA_main.pdf"), 
       p_pca_main, width = 10, height = 8, dpi = ANALYSIS_PARAMS$plot_dpi)
ggsave(file.path(PATHS$figures_qc, "PCA_main.png"), 
       p_pca_main, width = 10, height = 8, dpi = ANALYSIS_PARAMS$plot_dpi)

ggsave(file.path(PATHS$figures_qc, "PCA_by_group.pdf"), 
       p_pca_group, width = 10, height = 8, dpi = ANALYSIS_PARAMS$plot_dpi)

ggsave(file.path(PATHS$figures_qc, "PCA_PC3_PC4.pdf"), 
       p_pca_34, width = 10, height = 8, dpi = ANALYSIS_PARAMS$plot_dpi)

ggsave(file.path(PATHS$figures_qc, "scree_plot.pdf"), 
       p_scree, width = 8, height = 6, dpi = ANALYSIS_PARAMS$plot_dpi)

message("✓ PCA plots saved")

# =============================================================================
# 5. SAMPLE CORRELATION ANALYSIS
# =============================================================================

message("\n--- Sample Correlation Analysis ---\n")

# Calculate pairwise Pearson correlations
sample_cor <- cor(expr_matrix, method = "pearson")

message("Sample correlation matrix range: [", 
        round(min(sample_cor), 3), ", ", round(max(sample_cor), 3), "]")

# Check for samples with low correlation to others
mean_cor <- colMeans(sample_cor)
low_cor_threshold <- mean(mean_cor) - 2 * sd(mean_cor)

if (any(mean_cor < low_cor_threshold)) {
  low_cor_samples <- names(mean_cor)[mean_cor < low_cor_threshold]
  message("! Potential outlier samples (low correlation): ", 
          paste(low_cor_samples, collapse = ", "))
} else {
  message("✓ All samples show good correlation with others")
}

# Create annotation for heatmap
annotation_col <- sample_metadata[, c("cell_line", "treatment")]
rownames(annotation_col) <- sample_metadata$sample_id

# Define annotation colors
annotation_colors <- list(
  cell_line = ANALYSIS_PARAMS$colors$cell_line,
  treatment = ANALYSIS_PARAMS$colors$treatment
)

# Create correlation heatmap
pdf(file.path(PATHS$figures_qc, "sample_correlation_heatmap.pdf"), 
    width = 10, height = 8)
pheatmap(
  sample_cor,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks = seq(0.8, 1, length.out = 101),  # Focus on high correlation range
  main = "Sample Correlation Heatmap",
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 8,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black"
)
dev.off()

# Also save as PNG
png(file.path(PATHS$figures_qc, "sample_correlation_heatmap.png"), 
    width = 10, height = 8, units = "in", res = ANALYSIS_PARAMS$plot_dpi)
pheatmap(
  sample_cor,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks = seq(0.8, 1, length.out = 101),
  main = "Sample Correlation Heatmap",
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 8,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black"
)
dev.off()

message("✓ Correlation heatmap saved")

# =============================================================================
# 6. HIERARCHICAL CLUSTERING
# =============================================================================

message("\n--- Hierarchical Clustering ---\n")

# Perform hierarchical clustering on samples
sample_dist <- dist(t(expr_matrix), method = "euclidean")
sample_hclust <- hclust(sample_dist, method = "complete")

# Create dendrogram plot
pdf(file.path(PATHS$figures_qc, "sample_dendrogram.pdf"), width = 10, height = 6)
par(mar = c(8, 4, 4, 2))
plot(sample_hclust, 
     main = "Hierarchical Clustering of Samples",
     xlab = "",
     sub = "Complete linkage, Euclidean distance",
     cex = 0.8)

# Add colored rectangles for groups
# This helps visualize if samples cluster by biological groups
dev.off()

message("✓ Dendrogram saved")

# =============================================================================
# 7. EXPRESSION HEATMAP OF VARIABLE GENES
# =============================================================================

message("\n--- Variable Gene Expression Heatmap ---\n")

# Select top variable genes for visualization
gene_var <- apply(expr_matrix, 1, var)
top_var_genes <- names(sort(gene_var, decreasing = TRUE))[1:500]

# Subset expression matrix
expr_top_var <- expr_matrix[top_var_genes, ]

# Scale by row for visualization
expr_scaled <- t(scale(t(expr_top_var)))

# Create heatmap
pdf(file.path(PATHS$figures_qc, "heatmap_top_variable_genes.pdf"), 
    width = 10, height = 12)
pheatmap(
  expr_scaled,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-3, 3, length.out = 101),
  show_rownames = FALSE,  # Too many genes to show names
  main = "Top 500 Most Variable Genes",
  fontsize = 10,
  fontsize_col = 8
)
dev.off()

message("✓ Variable gene heatmap saved")

# =============================================================================
# 8. ASSESS BATCH EFFECTS AND TECHNICAL FACTORS
# =============================================================================

message("\n--- Assessing Technical Factors ---\n")

# Test association between PCs and experimental factors
pc_factor_assoc <- data.frame(
  Factor = c("Cell Line", "Treatment", "Group"),
  PC1_pvalue = NA,
  PC2_pvalue = NA,
  PC3_pvalue = NA
)

# PC1 associations
pc_factor_assoc$PC1_pvalue[1] <- summary(aov(PC1 ~ cell_line, data = pca_df))[[1]]$`Pr(>F)`[1]
pc_factor_assoc$PC1_pvalue[2] <- summary(aov(PC1 ~ treatment, data = pca_df))[[1]]$`Pr(>F)`[1]
pc_factor_assoc$PC1_pvalue[3] <- summary(aov(PC1 ~ group, data = pca_df))[[1]]$`Pr(>F)`[1]

# PC2 associations
pc_factor_assoc$PC2_pvalue[1] <- summary(aov(PC2 ~ cell_line, data = pca_df))[[1]]$`Pr(>F)`[1]
pc_factor_assoc$PC2_pvalue[2] <- summary(aov(PC2 ~ treatment, data = pca_df))[[1]]$`Pr(>F)`[1]
pc_factor_assoc$PC2_pvalue[3] <- summary(aov(PC2 ~ group, data = pca_df))[[1]]$`Pr(>F)`[1]

# PC3 associations  
pc_factor_assoc$PC3_pvalue[1] <- summary(aov(PC3 ~ cell_line, data = pca_df))[[1]]$`Pr(>F)`[1]
pc_factor_assoc$PC3_pvalue[2] <- summary(aov(PC3 ~ treatment, data = pca_df))[[1]]$`Pr(>F)`[1]
pc_factor_assoc$PC3_pvalue[3] <- summary(aov(PC3 ~ group, data = pca_df))[[1]]$`Pr(>F)`[1]

message("Association between PCs and experimental factors:")
print(pc_factor_assoc)

# Expected pattern:
# - PC1 should strongly associate with cell line (major source of variation)
# - Treatment effect may be captured in PC2 or PC3
# - If no association with biological factors, may indicate technical issues

# =============================================================================
# 9. OUTLIER DETECTION
# =============================================================================

message("\n--- Outlier Detection ---\n")

# Method 1: Based on sample distance from group centroid
outlier_scores <- data.frame(
  sample_id = sample_metadata$sample_id,
  group = sample_metadata$group,
  mean_correlation = mean_cor,
  distance_from_centroid = NA
)

# Calculate distance from group centroid
for (grp in levels(sample_metadata$group)) {
  grp_samples <- sample_metadata$sample_id[sample_metadata$group == grp]
  
  if (length(grp_samples) > 1) {
    grp_expr <- expr_matrix[, grp_samples]
    grp_centroid <- rowMeans(grp_expr)
    
    for (samp in grp_samples) {
      dist_to_centroid <- sqrt(sum((expr_matrix[, samp] - grp_centroid)^2))
      outlier_scores$distance_from_centroid[outlier_scores$sample_id == samp] <- dist_to_centroid
    }
  }
}

# Method 2: Based on median absolute deviation of expression
sample_mad <- apply(expr_matrix, 2, mad)
outlier_scores$mad <- sample_mad

# Flag potential outliers
# Use 2 SD threshold for each metric
outlier_scores$potential_outlier <- 
  (outlier_scores$mean_correlation < mean(outlier_scores$mean_correlation) - 2*sd(outlier_scores$mean_correlation)) |
  (outlier_scores$distance_from_centroid > mean(outlier_scores$distance_from_centroid, na.rm = TRUE) + 
     2*sd(outlier_scores$distance_from_centroid, na.rm = TRUE))

message("\nOutlier assessment:")
print(outlier_scores)

if (any(outlier_scores$potential_outlier)) {
  message("\n! Potential outliers detected: ", 
          paste(outlier_scores$sample_id[outlier_scores$potential_outlier], collapse = ", "))
  message("  Recommendation: Review carefully but proceed with caution given small sample size")
} else {
  message("\n✓ No obvious outliers detected")
}

# =============================================================================
# 10. CREATE QC SUMMARY REPORT
# =============================================================================

message("\n--- Creating QC Summary ---\n")

# Compile QC metrics
qc_summary <- list(
  # Data dimensions
  features = nrow(expr_matrix),
  samples = ncol(expr_matrix),
  
  # Distribution metrics
  expression_range = range(expr_matrix),
  sample_medians = sample_stats$median,
  median_cv = sd(sample_stats$median) / mean(sample_stats$median) * 100,
  
  # PCA results
  pca = pca_result,
  variance_explained = var_explained,
  pc_factor_associations = pc_factor_assoc,
  
  # Correlation
  sample_correlations = sample_cor,
  mean_within_group_cor = NA,  # Calculate below
  mean_between_group_cor = NA,
  
  # Outlier assessment
  outlier_scores = outlier_scores,
  potential_outliers = outlier_scores$sample_id[outlier_scores$potential_outlier]
)

# Calculate within-group and between-group correlations
within_cors <- c()
between_cors <- c()

for (i in 1:(ncol(sample_cor)-1)) {
  for (j in (i+1):ncol(sample_cor)) {
    samp_i <- colnames(sample_cor)[i]
    samp_j <- colnames(sample_cor)[j]
    
    grp_i <- sample_metadata$group[sample_metadata$sample_id == samp_i]
    grp_j <- sample_metadata$group[sample_metadata$sample_id == samp_j]
    
    if (grp_i == grp_j) {
      within_cors <- c(within_cors, sample_cor[i, j])
    } else {
      between_cors <- c(between_cors, sample_cor[i, j])
    }
  }
}

qc_summary$mean_within_group_cor <- mean(within_cors)
qc_summary$mean_between_group_cor <- mean(between_cors)

message("Within-group correlation (mean): ", round(qc_summary$mean_within_group_cor, 3))
message("Between-group correlation (mean): ", round(qc_summary$mean_between_group_cor, 3))

# Good: within-group should be higher than between-group
if (qc_summary$mean_within_group_cor > qc_summary$mean_between_group_cor) {
  message("✓ Samples are more similar within groups than between groups")
} else {
  message("! Unexpected: between-group correlation higher than within-group")
}

# Save QC output
qc_output_file <- file.path(PATHS$processed, "03_qc_output.rds")
saveRDS(qc_summary, qc_output_file)
message("\n✓ QC summary saved: ", qc_output_file)

# =============================================================================
# 11. GENERATE QC REPORT
# =============================================================================

qc_report <- paste0(
  "================================================================================\n",
  "QUALITY CONTROL REPORT\n",
  "================================================================================\n\n",
  "Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "Dataset: GSE74764\n\n",
  
  "DATA DIMENSIONS\n",
  "---------------\n",
  "Features (genes): ", nrow(expr_matrix), "\n",
  "Samples: ", ncol(expr_matrix), "\n\n",
  
  "EXPRESSION DISTRIBUTION\n",
  "-----------------------\n",
  "Range: [", round(min(expr_matrix), 2), ", ", round(max(expr_matrix), 2), "]\n",
  "Median CV across samples: ", round(qc_summary$median_cv, 2), "%\n",
  "Assessment: ", ifelse(qc_summary$median_cv < 10, "GOOD", "SUBOPTIMAL"), "\n\n",
  
  "PRINCIPAL COMPONENT ANALYSIS\n",
  "----------------------------\n",
  "PC1: ", round(var_explained[1], 1), "% variance\n",
  "PC2: ", round(var_explained[2], 1), "% variance\n",
  "PC3: ", round(var_explained[3], 1), "% variance\n",
  "Top 3 PCs explain: ", round(sum(var_explained[1:3]), 1), "% of variance\n\n",
  
  "PC-Factor Associations (p-values):\n",
  "  Cell Line: PC1=", format(pc_factor_assoc$PC1_pvalue[1], digits = 3), 
  ", PC2=", format(pc_factor_assoc$PC2_pvalue[1], digits = 3), "\n",
  "  Treatment: PC1=", format(pc_factor_assoc$PC1_pvalue[2], digits = 3), 
  ", PC2=", format(pc_factor_assoc$PC2_pvalue[2], digits = 3), "\n\n",
  
  "SAMPLE CORRELATIONS\n",
  "-------------------\n",
  "Mean within-group correlation: ", round(qc_summary$mean_within_group_cor, 3), "\n",
  "Mean between-group correlation: ", round(qc_summary$mean_between_group_cor, 3), "\n",
  "Assessment: ", ifelse(qc_summary$mean_within_group_cor > qc_summary$mean_between_group_cor, 
                         "GOOD - samples cluster by group", 
                         "CHECK - unexpected correlation pattern"), "\n\n",
  
  "OUTLIER ASSESSMENT\n",
  "------------------\n",
  "Potential outliers: ", ifelse(length(qc_summary$potential_outliers) == 0, 
                                  "None detected", 
                                  paste(qc_summary$potential_outliers, collapse = ", ")), "\n\n",
  
  "CONCLUSION\n",
  "----------\n",
  "Data quality appears ", 
  ifelse(qc_summary$median_cv < 10 && 
           qc_summary$mean_within_group_cor > qc_summary$mean_between_group_cor &&
           length(qc_summary$potential_outliers) == 0,
         "GOOD. Proceed with differential expression analysis.",
         "ACCEPTABLE with caveats. Review QC plots carefully."), "\n\n",
  
  "GENERATED FIGURES\n",
  "-----------------\n",
  "1. density_by_sample.pdf/png - Expression density plots\n",
  "2. boxplot_by_sample.pdf/png - Expression boxplots\n",
  "3. PCA_main.pdf/png - Principal component analysis\n",
  "4. sample_correlation_heatmap.pdf/png - Sample correlations\n",
  "5. sample_dendrogram.pdf - Hierarchical clustering\n",
  "6. heatmap_top_variable_genes.pdf - Variable gene heatmap\n",
  "7. scree_plot.pdf - PCA variance explained\n\n",
  
  "================================================================================\n"
)

# Save report
report_file <- file.path(PATHS$figures_qc, "QC_report.txt")
writeLines(qc_report, report_file)
message("✓ QC report saved: ", report_file)

# Print report
cat(qc_report)

# =============================================================================
# 12. CREATE COMBINED QC FIGURE (PUBLICATION-READY)
# =============================================================================

message("\n--- Creating Combined QC Figure ---\n")

# Create a publication-ready combined figure
combined_qc <- plot_grid(
  p_box_sample + theme(legend.position = "none"),
  p_pca_main + theme(legend.position = "bottom"),
  ncol = 1,
  labels = c("A", "B"),
  label_size = 14,
  rel_heights = c(0.4, 0.6)
)

ggsave(file.path(PATHS$figures_qc, "Figure_QC_combined.pdf"), 
       combined_qc, width = 12, height = 12, dpi = ANALYSIS_PARAMS$plot_dpi)
ggsave(file.path(PATHS$figures_qc, "Figure_QC_combined.png"), 
       combined_qc, width = 12, height = 12, dpi = ANALYSIS_PARAMS$plot_dpi)

message("✓ Combined QC figure saved")

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Quality control complete. Proceed to 04_differential_expression.R       ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

################################################################################
# END OF SCRIPT
################################################################################
