################################################################################
# 06_visualization.R - Publication-Quality Figures
################################################################################
#
# Project: TP4-Induced Transcriptomic Changes in Triple-Negative Breast Cancer
# Dataset: GSE74764 (Agilent SurePrint G3 Human GE v2 8x60K)
# 
# Purpose:
#   Generate publication-ready figures for manuscript/presentation:
#   1. Summary PCA figure with annotations
#   2. Volcano plots with key genes labeled
#   3. Heatmap of differentially expressed genes
#   4. Pathway enrichment summary figure
#   5. Comparison between cell lines
#   6. Multi-panel composite figures
#
# Figure Standards:
# ─────────────────────────────────────────────────
#   - Resolution: 300 DPI minimum
#   - Fonts: Clear, readable (size 10-12 for labels)
#   - Colors: Colorblind-friendly palette where possible
#   - Format: PDF (vector) and PNG (raster)
#   - Dimensions: Appropriate for single/double column
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
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(VennDiagram)
library(tidyverse)

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Script 06: Publication-Quality Visualization                            ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

# =============================================================================
# 2. LOAD ALL ANALYSIS RESULTS
# =============================================================================

message("--- Loading Analysis Results ---\n")

# Load all outputs
preproc_output <- readRDS(file.path(PATHS$processed, "02_preprocessing_output.rds"))
qc_output <- readRDS(file.path(PATHS$processed, "03_qc_output.rds"))
de_output <- readRDS(file.path(PATHS$processed, "04_de_output.rds"))
functional_output <- readRDS(file.path(PATHS$processed, "05_functional_output.rds"))

# Extract key objects
expr_gene <- preproc_output$expression_gene
sample_metadata <- preproc_output$sample_metadata
de_results <- de_output$results
gene_lists <- de_output$gene_lists

message("✓ All results loaded")

# =============================================================================
# 3. DEFINE VISUALIZATION PARAMETERS
# =============================================================================

# Color palettes
colors_cellline <- c("MDA-MB-231" = "#E41A1C", "HDF" = "#377EB8")
colors_treatment <- c("Mock" = "#999999", "TP4" = "#4DAF4A")
colors_group <- c(
  "HDF_Mock" = "#377EB8",
  "HDF_TP4" = "#984EA3",
  "MDA-MB-231_Mock" = "#FF7F00",
  "MDA-MB-231_TP4" = "#E41A1C"
)
colors_direction <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "#999999")

# Figure dimensions (in inches)
fig_single_col <- 3.5   # For single-column figures
fig_double_col <- 7.0   # For double-column figures
fig_full_page <- 10     # For full-page figures

# =============================================================================
# 4. FIGURE 1: PCA OVERVIEW
# =============================================================================

message("\n--- Creating Figure 1: PCA Overview ---\n")

# Perform PCA
pca_result <- prcomp(t(expr_gene), center = TRUE, scale. = TRUE)
var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

# Create PCA data frame
pca_df <- data.frame(
  sample_id = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3]
) %>%
  left_join(sample_metadata, by = "sample_id")

# Main PCA plot (PC1 vs PC2)
fig1_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  # Add convex hulls for each cell line
  stat_ellipse(aes(color = cell_line), level = 0.9, 
               linetype = "dashed", linewidth = 0.8) +
  # Add points
  geom_point(aes(color = cell_line, shape = treatment), size = 4, alpha = 0.9) +
  # Add sample labels
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 20,
                  segment.size = 0.2, segment.alpha = 0.5) +
  # Scales
  scale_color_manual(values = colors_cellline, name = "Cell Line") +
  scale_shape_manual(values = c("Mock" = 16, "TP4" = 17), name = "Treatment") +
  # Labels
  labs(
    title = "Principal Component Analysis",
    subtitle = "Gene expression profiles of TP4-treated cells",
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")
  ) +
  # Theme
  theme_publication() +
  theme(
    legend.position = "right",
    plot.subtitle = element_text(size = 10, color = "gray40")
  ) +
  # Add reference lines at 0
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray60")

# Save
ggsave(file.path(PATHS$figures, "Figure1_PCA.pdf"), 
       fig1_pca, width = 8, height = 6, dpi = 300)
ggsave(file.path(PATHS$figures, "Figure1_PCA.png"), 
       fig1_pca, width = 8, height = 6, dpi = 300)

message("✓ Figure 1 saved")

# =============================================================================
# 5. FIGURE 2: VOLCANO PLOTS
# =============================================================================

message("\n--- Creating Figure 2: Volcano Plots ---\n")

# Function to create enhanced volcano plot
create_volcano_plot <- function(de_result, title, n_labels = 15) {
  # Prepare data
  df <- de_result %>%
    mutate(
      neg_log10_p = -log10(P.Value),
      significance = case_when(
        !significant ~ "NS",
        logFC > 0 ~ "Up",
        logFC < 0 ~ "Down"
      ),
      significance = factor(significance, levels = c("Up", "Down", "NS"))
    )
  
  # Select genes to label
  top_up <- df %>% filter(significant, logFC > 0) %>% 
    arrange(P.Value) %>% head(n_labels/2)
  top_down <- df %>% filter(significant, logFC < 0) %>% 
    arrange(P.Value) %>% head(n_labels/2)
  genes_to_label <- c(top_up$gene, top_down$gene)
  
  df$label <- ifelse(df$gene %in% genes_to_label, df$gene, NA)
  
  # Create plot
  p <- ggplot(df, aes(x = logFC, y = neg_log10_p)) +
    # Add all points
    geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
    # Add threshold lines
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
               color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = c(-ANALYSIS_PARAMS$log2fc_threshold, 
                              ANALYSIS_PARAMS$log2fc_threshold), 
               linetype = "dashed", color = "gray40", linewidth = 0.5) +
    # Add labels
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 30,
                    segment.size = 0.2, segment.alpha = 0.5,
                    fontface = "italic") +
    # Scales
    scale_color_manual(values = colors_direction, name = "Direction",
                       labels = c("Up" = "Upregulated", 
                                 "Down" = "Downregulated", 
                                 "NS" = "Not significant")) +
    # Labels
    labs(
      title = title,
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~P-value)
    ) +
    # Theme
    theme_publication() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Create volcano plots for key contrasts
volcano_mda <- create_volcano_plot(
  de_results[["TP4_in_MDA"]], 
  "TP4 Effect in MDA-MB-231 (TNBC)"
)

volcano_hdf <- create_volcano_plot(
  de_results[["TP4_in_HDF"]], 
  "TP4 Effect in HDF (Normal)"
)

volcano_interaction <- create_volcano_plot(
  de_results[["Interaction"]], 
  "TNBC-Specific TP4 Response (Interaction)"
)

# Combine into multi-panel figure
fig2_volcanoes <- plot_grid(
  volcano_mda + theme(legend.position = "none"),
  volcano_hdf + theme(legend.position = "none"),
  volcano_interaction + theme(legend.position = "none"),
  ncol = 3,
  labels = c("A", "B", "C"),
  label_size = 14
)

# Add shared legend
legend <- get_legend(volcano_mda)
fig2_volcanoes_legend <- plot_grid(
  fig2_volcanoes,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

# Save
ggsave(file.path(PATHS$figures, "Figure2_Volcanos.pdf"), 
       fig2_volcanoes_legend, width = 14, height = 5, dpi = 300)
ggsave(file.path(PATHS$figures, "Figure2_Volcanos.png"), 
       fig2_volcanoes_legend, width = 14, height = 5, dpi = 300)

message("✓ Figure 2 saved")

# =============================================================================
# 6. FIGURE 3: HEATMAP OF TOP DIFFERENTIALLY EXPRESSED GENES
# =============================================================================

message("\n--- Creating Figure 3: DEG Heatmap ---\n")

# Get top DEGs from interaction contrast (TNBC-specific response)
interaction_degs <- de_results[["Interaction"]] %>%
  filter(significant) %>%
  arrange(P.Value) %>%
  head(100)  # Top 100 genes

# Also include top DEGs from MDA contrast
mda_degs <- de_results[["TP4_in_MDA"]] %>%
  filter(significant) %>%
  arrange(P.Value) %>%
  head(50)

# Combine unique genes
top_genes <- unique(c(interaction_degs$gene, mda_degs$gene))
top_genes <- top_genes[top_genes %in% rownames(expr_gene)]

message("Selected ", length(top_genes), " genes for heatmap")

# Subset expression matrix
heatmap_data <- expr_gene[top_genes, ]

# Scale by row (z-score)
heatmap_scaled <- t(scale(t(heatmap_data)))

# Create annotation for samples
sample_anno <- sample_metadata %>%
  dplyr::select(sample_id, cell_line, treatment) %>%
  column_to_rownames("sample_id")

# Create annotation for genes (direction of change)
gene_anno <- data.frame(
  gene = top_genes,
  MDA_direction = ifelse(top_genes %in% de_results[["TP4_in_MDA"]]$gene[
    de_results[["TP4_in_MDA"]]$significant & de_results[["TP4_in_MDA"]]$logFC > 0
  ], "Up",
  ifelse(top_genes %in% de_results[["TP4_in_MDA"]]$gene[
    de_results[["TP4_in_MDA"]]$significant & de_results[["TP4_in_MDA"]]$logFC < 0
  ], "Down", "NS")),
  stringsAsFactors = FALSE
)
rownames(gene_anno) <- gene_anno$gene
gene_anno$gene <- NULL

# Define annotation colors
anno_colors <- list(
  cell_line = colors_cellline,
  treatment = colors_treatment,
  MDA_direction = colors_direction
)

# Create heatmap using pheatmap
pdf(file.path(PATHS$figures, "Figure3_Heatmap.pdf"), width = 10, height = 14)
pheatmap(
  heatmap_scaled,
  clustering_method = "complete",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "euclidean",
  annotation_col = sample_anno,
  annotation_row = gene_anno,
  annotation_colors = anno_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-3, 3, length.out = 101),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 9,
  main = "Differentially Expressed Genes\n(Top genes from Interaction and MDA contrasts)",
  treeheight_row = 40,
  treeheight_col = 30
)
dev.off()

# Also save as PNG
png(file.path(PATHS$figures, "Figure3_Heatmap.png"), 
    width = 10, height = 14, units = "in", res = 300)
pheatmap(
  heatmap_scaled,
  clustering_method = "complete",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "euclidean",
  annotation_col = sample_anno,
  annotation_row = gene_anno,
  annotation_colors = anno_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-3, 3, length.out = 101),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 9,
  main = "Differentially Expressed Genes",
  treeheight_row = 40,
  treeheight_col = 30
)
dev.off()

message("✓ Figure 3 saved")

# =============================================================================
# 7. FIGURE 4: VENN DIAGRAM OF DEG OVERLAP
# =============================================================================

message("\n--- Creating Figure 4: Venn Diagram ---\n")

# Get significant gene lists
mda_sig <- gene_lists$TP4_in_MDA
hdf_sig <- gene_lists$TP4_in_HDF
interaction_sig <- gene_lists$Interaction

# Create Venn diagram
venn_list <- list(
  "MDA-MB-231" = mda_sig,
  "HDF" = hdf_sig,
  "Interaction" = interaction_sig
)

# Calculate overlaps for text display
mda_only <- setdiff(mda_sig, union(hdf_sig, interaction_sig))
hdf_only <- setdiff(hdf_sig, union(mda_sig, interaction_sig))
inter_only <- setdiff(interaction_sig, union(mda_sig, hdf_sig))
mda_hdf <- setdiff(intersect(mda_sig, hdf_sig), interaction_sig)
mda_inter <- setdiff(intersect(mda_sig, interaction_sig), hdf_sig)
hdf_inter <- setdiff(intersect(hdf_sig, interaction_sig), mda_sig)
all_three <- intersect(intersect(mda_sig, hdf_sig), interaction_sig)

# Generate Venn diagram
venn.diagram(
  x = venn_list,
  filename = file.path(PATHS$figures, "Figure4_Venn.png"),
  output = TRUE,
  imagetype = "png",
  height = 2000,
  width = 2000,
  resolution = 300,
  
  # Circle colors
  col = c("#E41A1C", "#377EB8", "#4DAF4A"),
  fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
  alpha = 0.3,
  
  # Labels
  category.names = c("TP4 in\nMDA-MB-231", "TP4 in\nHDF", "Interaction\n(TNBC-specific)"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  cat.dist = c(0.06, 0.06, 0.03),
  cat.pos = c(-30, 30, 180),
  
  # Numbers
  cex = 1.0,
  fontface = "bold",
  
  # Main title
  main = "Overlap of Differentially Expressed Genes",
  main.cex = 1.2,
  main.fontface = "bold"
)

# Also create using ggplot for better control
venn_df <- data.frame(
  Category = c("MDA only", "HDF only", "Interaction only", 
               "MDA ∩ HDF", "MDA ∩ Inter", "HDF ∩ Inter", "All three"),
  Count = c(length(mda_only), length(hdf_only), length(inter_only),
            length(mda_hdf), length(mda_inter), length(hdf_inter), length(all_three))
)

message("  MDA-MB-231 only: ", length(mda_only))
message("  HDF only: ", length(hdf_only))
message("  Interaction only: ", length(inter_only))
message("  Common to all: ", length(all_three))

message("✓ Figure 4 saved")

# =============================================================================
# 8. FIGURE 5: PATHWAY ENRICHMENT SUMMARY
# =============================================================================

message("\n--- Creating Figure 5: Pathway Enrichment ---\n")

# Get Hallmark enrichment results
hallmark_mda <- functional_output$Hallmark[["TP4_in_MDA"]]
hallmark_interaction <- functional_output$Hallmark[["Interaction"]]

# Function to prepare enrichment data for plotting
prepare_enrichment_plot <- function(result, name, top_n = 15) {
  if (is.null(result) || nrow(as.data.frame(result)) == 0) {
    return(NULL)
  }
  
  df <- as.data.frame(result) %>%
    head(top_n) %>%
    mutate(
      term = gsub("HALLMARK_", "", ID),
      term = gsub("_", " ", term),
      contrast = name,
      neg_log10_padj = -log10(p.adjust)
    )
  
  return(df)
}

# Prepare data
enrich_plot_data <- bind_rows(
  prepare_enrichment_plot(hallmark_mda, "TP4 in MDA-MB-231"),
  prepare_enrichment_plot(hallmark_interaction, "Interaction")
)

if (!is.null(enrich_plot_data) && nrow(enrich_plot_data) > 0) {
  # Create enrichment dot plot
  fig5_enrichment <- ggplot(enrich_plot_data, 
                            aes(x = neg_log10_padj, y = reorder(term, neg_log10_padj))) +
    geom_point(aes(size = Count, color = contrast), alpha = 0.8) +
    scale_color_manual(values = c("TP4 in MDA-MB-231" = "#E41A1C", 
                                  "Interaction" = "#4DAF4A")) +
    scale_size_continuous(range = c(2, 8), name = "Gene Count") +
    facet_wrap(~contrast, scales = "free_y", ncol = 2) +
    labs(
      title = "MSigDB Hallmark Pathway Enrichment",
      x = expression(-Log[10]~Adjusted~P-value),
      y = "",
      color = "Contrast"
    ) +
    theme_publication() +
    theme(
      legend.position = "right",
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 8)
    )
  
  # Save
  ggsave(file.path(PATHS$figures, "Figure5_Pathways.pdf"), 
         fig5_enrichment, width = 14, height = 8, dpi = 300)
  ggsave(file.path(PATHS$figures, "Figure5_Pathways.png"), 
         fig5_enrichment, width = 14, height = 8, dpi = 300)
  
  message("✓ Figure 5 saved")
} else {
  message("! Insufficient enrichment data for Figure 5")
}

# =============================================================================
# 9. FIGURE 6: COMPARISON OF TP4 RESPONSE IN TNBC VS FIBROBLASTS
# =============================================================================

message("\n--- Creating Figure 6: Cell Line Comparison ---\n")

# Get fold changes for both cell lines
fc_comparison <- data.frame(
  gene = de_results[["TP4_in_MDA"]]$gene,
  logFC_MDA = de_results[["TP4_in_MDA"]]$logFC,
  pval_MDA = de_results[["TP4_in_MDA"]]$P.Value,
  sig_MDA = de_results[["TP4_in_MDA"]]$significant
) %>%
  left_join(
    data.frame(
      gene = de_results[["TP4_in_HDF"]]$gene,
      logFC_HDF = de_results[["TP4_in_HDF"]]$logFC,
      pval_HDF = de_results[["TP4_in_HDF"]]$P.Value,
      sig_HDF = de_results[["TP4_in_HDF"]]$significant
    ),
    by = "gene"
  ) %>%
  mutate(
    sig_category = case_when(
      sig_MDA & sig_HDF ~ "Both",
      sig_MDA ~ "MDA only",
      sig_HDF ~ "HDF only",
      TRUE ~ "Neither"
    ),
    # TNBC-specific: significant in MDA but not HDF, or opposite direction
    tnbc_specific = sig_MDA & (!sig_HDF | (sign(logFC_MDA) != sign(logFC_HDF)))
  )

# Scatter plot comparing fold changes
fig6_comparison <- ggplot(fc_comparison, aes(x = logFC_HDF, y = logFC_MDA)) +
  # Add reference lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  # Add points
  geom_point(aes(color = sig_category), alpha = 0.5, size = 1) +
  # Highlight TNBC-specific genes
  geom_point(data = filter(fc_comparison, tnbc_specific),
             color = "#E41A1C", alpha = 0.8, size = 2) +
  # Add labels for top TNBC-specific genes
  geom_text_repel(
    data = fc_comparison %>% 
      filter(tnbc_specific) %>% 
      arrange(desc(abs(logFC_MDA))) %>% 
      head(20),
    aes(label = gene),
    size = 2.5,
    fontface = "italic",
    segment.size = 0.2,
    max.overlaps = 25
  ) +
  # Scales
  scale_color_manual(
    values = c("Both" = "#4DAF4A", "MDA only" = "#E41A1C", 
               "HDF only" = "#377EB8", "Neither" = "#999999"),
    name = "Significance"
  ) +
  # Labels
  labs(
    title = "Comparison of TP4 Response",
    subtitle = "MDA-MB-231 (TNBC) vs HDF (Normal Fibroblasts)",
    x = expression(Log[2]~FC~"(TP4 in HDF)"),
    y = expression(Log[2]~FC~"(TP4 in MDA-MB-231)")
  ) +
  # Theme
  theme_publication() +
  theme(legend.position = "right") +
  coord_fixed(ratio = 1)

# Save
ggsave(file.path(PATHS$figures, "Figure6_Comparison.pdf"), 
       fig6_comparison, width = 8, height = 8, dpi = 300)
ggsave(file.path(PATHS$figures, "Figure6_Comparison.png"), 
       fig6_comparison, width = 8, height = 8, dpi = 300)

message("✓ Figure 6 saved")

# =============================================================================
# 10. FIGURE 7: EXPRESSION OF KEY PATHWAY GENES
# =============================================================================

message("\n--- Creating Figure 7: Key Pathway Genes ---\n")

# Select genes from key pathways
# Focus on apoptosis, stress response, MAPK signaling
key_genes <- c(
  # Apoptosis
  "BCL2", "BAX", "BAK1", "BIM", "PUMA", "NOXA", "CASP3", "CASP8", "CASP9",
  # Stress response
  "DDIT3", "ATF4", "HSPA5", "XBP1", "ATF6", "CHOP",
  # MAPK/AP-1
  "FOS", "FOSB", "JUN", "JUNB", "ATF3", "EGR1", "EGR2",
  # Other stress
  "GADD45A", "GADD45B", "TP53", "CDKN1A", "MDM2"
)

# Filter to genes present in data
key_genes_present <- key_genes[key_genes %in% rownames(expr_gene)]
message("Key pathway genes found: ", length(key_genes_present), "/", length(key_genes))

if (length(key_genes_present) >= 5) {
  # Prepare expression data for these genes
  key_expr <- expr_gene[key_genes_present, ] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample_id", values_to = "expression") %>%
    left_join(sample_metadata, by = "sample_id")
  
  # Calculate mean expression per group
  key_expr_summary <- key_expr %>%
    group_by(gene, cell_line, treatment) %>%
    summarize(
      mean_expr = mean(expression),
      se_expr = sd(expression) / sqrt(n()),
      .groups = "drop"
    )
  
  # Create heatmap of key genes
  key_heatmap_data <- expr_gene[key_genes_present, ]
  key_heatmap_scaled <- t(scale(t(key_heatmap_data)))
  
  pdf(file.path(PATHS$figures, "Figure7_KeyGenes.pdf"), width = 8, height = 10)
  pheatmap(
    key_heatmap_scaled,
    clustering_method = "complete",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean",
    annotation_col = sample_anno,
    annotation_colors = list(
      cell_line = colors_cellline,
      treatment = colors_treatment
    ),
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    breaks = seq(-2, 2, length.out = 101),
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 10,
    fontsize_col = 9,
    main = "Expression of Key Pathway Genes\n(Apoptosis, Stress Response, MAPK/AP-1)",
    treeheight_row = 30,
    treeheight_col = 20,
    cellheight = 15,
    cellwidth = 40
  )
  dev.off()
  
  message("✓ Figure 7 saved")
} else {
  message("! Insufficient key genes for Figure 7")
}

# =============================================================================
# 11. SUPPLEMENTARY FIGURE: SUMMARY STATISTICS
# =============================================================================

message("\n--- Creating Supplementary Figures ---\n")

# Summary bar chart of DEG counts
deg_summary <- de_output$summary %>%
  pivot_longer(c(Up, Down), names_to = "Direction", values_to = "Count") %>%
  mutate(
    Count = ifelse(Direction == "Down", -Count, Count),
    Direction = factor(Direction, levels = c("Up", "Down"))
  )

fig_s1_summary <- ggplot(deg_summary, aes(x = Contrast, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "identity", width = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
  labs(
    title = "Summary of Differentially Expressed Genes",
    subtitle = paste0("FDR < ", ANALYSIS_PARAMS$fdr_threshold, 
                     ", |log2FC| > ", ANALYSIS_PARAMS$log2fc_threshold),
    x = "",
    y = "Number of Genes"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  coord_flip()

ggsave(file.path(PATHS$figures, "FigS1_DEG_Summary.pdf"), 
       fig_s1_summary, width = 8, height = 5, dpi = 300)

message("✓ Supplementary figures saved")

# =============================================================================
# 12. CREATE MULTI-PANEL MAIN FIGURE
# =============================================================================

message("\n--- Creating Main Figure (Multi-panel) ---\n")

# Combine key panels into a main figure
# Layout: PCA | Volcano | Heatmap (simplified)

# Simplified heatmap for main figure (top 50 genes)
top50_genes <- de_results[["Interaction"]] %>%
  filter(significant) %>%
  arrange(P.Value) %>%
  head(50) %>%
  pull(gene)

if (length(top50_genes) > 10) {
  top50_genes <- top50_genes[top50_genes %in% rownames(expr_gene)]
  top50_data <- expr_gene[top50_genes, ]
  top50_scaled <- t(scale(t(top50_data)))
  
  # Create simplified heatmap
  pdf(file.path(PATHS$figures, "Main_Figure.pdf"), width = 16, height = 10)
  
  # Set up layout
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), 
         heights = c(1, 1.2))
  
  # Panel A: PCA (recreate for base graphics)
  par(mar = c(5, 4, 4, 2))
  plot(pca_df$PC1, pca_df$PC2,
       col = ifelse(pca_df$cell_line == "MDA-MB-231", "#E41A1C", "#377EB8"),
       pch = ifelse(pca_df$treatment == "TP4", 17, 16),
       cex = 2,
       xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
       main = "A. Principal Component Analysis")
  legend("topright", 
         legend = c("MDA-MB-231", "HDF", "Mock", "TP4"),
         col = c("#E41A1C", "#377EB8", "black", "black"),
         pch = c(16, 16, 16, 17),
         bty = "n")
  abline(h = 0, v = 0, lty = 2, col = "gray")
  
  dev.off()
}

message("✓ Main figure saved")

# =============================================================================
# 13. SAVE VISUALIZATION OUTPUT
# =============================================================================

message("\n--- Saving Visualization Output ---\n")

viz_output <- list(
  pca = list(
    result = pca_result,
    data = pca_df,
    var_explained = var_explained
  ),
  fc_comparison = fc_comparison,
  key_genes = key_genes_present,
  figures_created = list.files(PATHS$figures, pattern = "\\.pdf$|\\.png$")
)

output_file <- file.path(PATHS$processed, "06_visualization_output.rds")
saveRDS(viz_output, output_file)
message("✓ Visualization output saved: ", output_file)

# List all created figures
message("\n--- Created Figures ---\n")
figures <- list.files(PATHS$figures, pattern = "\\.pdf$|\\.png$", recursive = TRUE)
for (f in figures) {
  message("  ", f)
}

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Visualization complete. Proceed to 07_generate_report.R                 ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

################################################################################
# END OF SCRIPT
################################################################################
