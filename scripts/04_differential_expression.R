################################################################################
# 04_differential_expression.R - Differential Expression Analysis with limma
################################################################################
#
# Project: TP4-Induced Transcriptomic Changes in Triple-Negative Breast Cancer
# Dataset: GSE74764 (Agilent SurePrint G3 Human GE v2 8x60K)
# 
# Purpose:
#   Perform rigorous differential expression analysis using limma:
#   1. Build appropriate design matrix for 2x2 factorial design
#   2. Define biologically meaningful contrasts
#   3. Fit linear models with empirical Bayes moderation
#   4. Extract and annotate differentially expressed genes
#   5. Generate diagnostic plots
#
# Statistical Model:
# ─────────────────────────────────────────────────
#   We use a CELL-MEANS MODEL (group model) for clarity:
#   
#   Expression ~ 0 + Group
#   
#   Where Group = {HDF_Mock, HDF_TP4, MDA_Mock, MDA_TP4}
#   
#   This parameterization allows straightforward definition of contrasts:
#   
#   CONTRAST A: TP4 effect in MDA-MB-231 (TNBC)
#              = MDA_TP4 - MDA_Mock
#   
#   CONTRAST B: TP4 effect in HDF (normal fibroblasts)
#              = HDF_TP4 - HDF_Mock
#   
#   CONTRAST C: Baseline difference (cancer vs normal)
#              = MDA_Mock - HDF_Mock
#   
#   CONTRAST D: INTERACTION (cancer-specific TP4 response)
#              = (MDA_TP4 - MDA_Mock) - (HDF_TP4 - HDF_Mock)
#              
#   The INTERACTION contrast (D) is the most scientifically important
#   because it identifies genes that respond differently to TP4 in
#   TNBC cells compared to normal fibroblasts - these are candidates
#   for explaining TP4's selective anticancer activity.
#
# Rationale for Statistical Approach:
# ─────────────────────────────────────────────────
#   1. EMPIRICAL BAYES (eBayes):
#      With only n=2 per group, gene-specific variance estimates are
#      unreliable. limma's empirical Bayes approach borrows information
#      across genes to stabilize variance estimates, improving power
#      while controlling false positives.
#   
#   2. FDR THRESHOLD (0.05):
#      Standard threshold for discovery; given small sample size,
#      we emphasize effect sizes alongside statistical significance.
#   
#   3. FOLD-CHANGE THRESHOLD (|log2FC| > 0.585 = 1.5-fold):
#      Less stringent than the typical 2-fold cutoff to capture
#      biologically meaningful responses despite limited replication.
#
# Author: Computational Biology Pipeline
# Date: 2026-01-19
#
################################################################################

# =============================================================================
# 1. SETUP (Windows/RStudio Compatible)
# =============================================================================

# Find and source the setup script (works from any working directory)
local({
  paths_to_try <- c(
    "scripts/00_setup.R",
    "00_setup.R",
    "../scripts/00_setup.R"
  )
  if (requireNamespace("here", quietly = TRUE)) {
    paths_to_try <- c(here::here("scripts", "00_setup.R"), paths_to_try)
  }
  for (p in paths_to_try) {
    if (file.exists(p)) { source(p); return() }
  }
  stop("Cannot find 00_setup.R. Please set working directory to project root.")
})

# Load required packages
library(limma)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Script 04: Differential Expression Analysis                             ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

# =============================================================================
# 2. LOAD PREPROCESSED DATA
# =============================================================================

message("--- Loading Preprocessed Data ---\n")

# Load preprocessing output
preproc_output <- readRDS(file.path(PATHS$processed, "02_preprocessing_output.rds"))

expr_gene <- preproc_output$expression_gene
expr_probe <- preproc_output$expression_probe
sample_metadata <- preproc_output$sample_metadata

# Use gene-level expression for main analysis
if (!is.null(expr_gene)) {
  expr_matrix <- expr_gene
  analysis_level <- "gene"
  message("Using gene-level expression matrix: ", nrow(expr_matrix), " genes")
} else {
  expr_matrix <- expr_probe
  analysis_level <- "probe"
  message("Using probe-level expression matrix: ", nrow(expr_matrix), " probes")
}

message("Samples: ", ncol(expr_matrix))

# Verify sample order matches
stopifnot(all(colnames(expr_matrix) == sample_metadata$sample_id))
message("✓ Sample order verified")

# =============================================================================
# 3. DEFINE EXPERIMENTAL DESIGN
# =============================================================================

message("\n--- Defining Experimental Design ---\n")

# Display the experimental design
message("Experimental design (2x2 factorial):")
print(table(sample_metadata$cell_line, sample_metadata$treatment))

# Ensure factors are properly coded
sample_metadata$group <- factor(
  sample_metadata$group,
  levels = c("HDF_Mock", "HDF_TP4", "MDA-MB-231_Mock", "MDA-MB-231_TP4")
)

# Check for balanced design
group_counts <- table(sample_metadata$group)
message("\nSamples per group:")
print(group_counts)

if (any(group_counts < 2)) {
  warning("Some groups have fewer than 2 replicates!")
}

# =============================================================================
# 4. BUILD DESIGN MATRIX
# =============================================================================

message("\n--- Building Design Matrix ---\n")

# CELL-MEANS MODEL (no intercept)
# This parameterization estimates the mean of each group directly,
# making contrast specification more intuitive

design <- model.matrix(~ 0 + group, data = sample_metadata)

# Clean up column names (remove "group" prefix)
colnames(design) <- gsub("group", "", colnames(design))
# Handle special characters in column names
colnames(design) <- make.names(colnames(design))

message("Design matrix:")
print(design)

# Verify design matrix properties
message("\nDesign matrix dimensions: ", nrow(design), " samples × ", ncol(design), " coefficients")
message("Column names: ", paste(colnames(design), collapse = ", "))

# Check for rank deficiency
if (qr(design)$rank < ncol(design)) {
  warning("Design matrix is rank-deficient!")
} else {
  message("✓ Design matrix has full rank")
}

# =============================================================================
# 5. DEFINE CONTRASTS
# =============================================================================

message("\n--- Defining Contrasts ---\n")

# Get actual column names from design matrix
col_names <- colnames(design)
message("Design matrix columns: ", paste(col_names, collapse = ", "))

# Build contrast matrix
# Note: Column names may have dots instead of hyphens due to make.names()

contrast_matrix <- makeContrasts(
  # CONTRAST A: TP4 effect in MDA-MB-231 (TNBC)
  # Genes changed by TP4 treatment in cancer cells
  TP4_in_MDA = MDA.MB.231_TP4 - MDA.MB.231_Mock,
  
  # CONTRAST B: TP4 effect in HDF (normal fibroblasts)
  # Genes changed by TP4 treatment in normal cells
  TP4_in_HDF = HDF_TP4 - HDF_Mock,
  
  # CONTRAST C: Baseline cancer vs normal (without treatment)
  # Inherent differences between cancer and normal cells
  MDA_vs_HDF_baseline = MDA.MB.231_Mock - HDF_Mock,
  
  # CONTRAST D: INTERACTION - TNBC-specific TP4 response
  # This is the KEY contrast: identifies genes that respond
  # differently to TP4 in cancer cells vs normal cells
  # Positive: more upregulated in MDA by TP4 than in HDF
  # Negative: more downregulated in MDA by TP4 than in HDF
  Interaction = (MDA.MB.231_TP4 - MDA.MB.231_Mock) - (HDF_TP4 - HDF_Mock),
  
  levels = design
)

message("Contrast matrix:")
print(contrast_matrix)

message("\nContrast interpretations:")
message("  TP4_in_MDA:          Genes changed by TP4 in TNBC cells")
message("  TP4_in_HDF:          Genes changed by TP4 in normal fibroblasts")
message("  MDA_vs_HDF_baseline: Baseline cancer vs normal differences")
message("  Interaction:         TNBC-SPECIFIC TP4 response (key contrast)")

# =============================================================================
# 6. FIT LINEAR MODEL
# =============================================================================

message("\n--- Fitting Linear Model ---\n")

# Step 1: Fit the linear model
fit <- lmFit(expr_matrix, design)
message("✓ Linear model fitted")

# Step 2: Apply contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
message("✓ Contrasts applied")

# Step 3: Apply empirical Bayes moderation
# This is CRITICAL for small sample sizes
# It borrows information across genes to stabilize variance estimates
fit2 <- eBayes(fit2)
message("✓ Empirical Bayes moderation applied")

# Display summary of the fit
message("\nModel summary:")
message("  Genes analyzed: ", nrow(fit2$coefficients))
message("  Contrasts: ", ncol(fit2$coefficients))
message("  Prior degrees of freedom: ", round(fit2$df.prior, 2))
message("  Prior variance: ", round(fit2$s2.prior, 4))

# The prior df indicates how much information is borrowed
# Higher values = more borrowing = more stable estimates
if (fit2$df.prior > 3) {
  message("  ✓ Substantial variance shrinkage achieved (df.prior > 3)")
} else {
  message("  ! Limited variance shrinkage - interpret with caution")
}

# =============================================================================
# 7. EXTRACT DIFFERENTIAL EXPRESSION RESULTS
# =============================================================================

message("\n--- Extracting DE Results ---\n")

# Define thresholds
fdr_threshold <- ANALYSIS_PARAMS$fdr_threshold
log2fc_threshold <- ANALYSIS_PARAMS$log2fc_threshold

message("Significance thresholds:")
message("  FDR < ", fdr_threshold)
message("  |log2FC| > ", log2fc_threshold, " (", round(2^log2fc_threshold, 2), "-fold)")

# Extract results for each contrast
contrast_names <- colnames(contrast_matrix)
de_results <- list()
de_summary <- data.frame(
  Contrast = contrast_names,
  Total_Tested = nrow(fit2$coefficients),
  Sig_FDR = NA,
  Sig_FDR_FC = NA,
  Up = NA,
  Down = NA
)

for (i in seq_along(contrast_names)) {
  contrast <- contrast_names[i]
  message("\n--- ", contrast, " ---")
  
  # Extract topTable with all genes
  tt <- topTable(
    fit2,
    coef = contrast,
    number = Inf,           # All genes
    adjust.method = "BH",   # Benjamini-Hochberg FDR
    sort.by = "P"
  )
  
  # Add gene names if rownames are genes
  tt$gene <- rownames(tt)
  
  # Add direction and significance
  tt$direction <- ifelse(tt$logFC > 0, "Up", "Down")
  tt$significant_fdr <- tt$adj.P.Val < fdr_threshold
  tt$significant_fc <- abs(tt$logFC) > log2fc_threshold
  tt$significant <- tt$significant_fdr & tt$significant_fc
  
  # Calculate fold change from log2FC
  tt$fold_change <- 2^tt$logFC
  tt$fold_change_direction <- ifelse(tt$logFC > 0, tt$fold_change, -1/tt$fold_change)
  
  # Count significant genes
  n_sig_fdr <- sum(tt$significant_fdr)
  n_sig_both <- sum(tt$significant)
  n_up <- sum(tt$significant & tt$direction == "Up")
  n_down <- sum(tt$significant & tt$direction == "Down")
  
  message("  FDR < ", fdr_threshold, ": ", n_sig_fdr, " genes")
  message("  FDR + FC threshold: ", n_sig_both, " genes (", n_up, " up, ", n_down, " down)")
  
  # Store results
  de_results[[contrast]] <- tt
  de_summary$Sig_FDR[i] <- n_sig_fdr
  de_summary$Sig_FDR_FC[i] <- n_sig_both
  de_summary$Up[i] <- n_up
  de_summary$Down[i] <- n_down
  
  # Show top genes
  message("  Top 5 genes by significance:")
  print(head(tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], 5))
}

# Display summary
message("\n--- DE Summary ---\n")
print(de_summary)

# =============================================================================
# 8. ANNOTATE RESULTS WITH BIOLOGICAL CONTEXT
# =============================================================================

message("\n--- Annotating Results ---\n")

# Try to add additional annotations using org.Hs.eg.db
# Note: This may not work if gene symbols are non-standard

library(org.Hs.eg.db)
library(AnnotationDbi)

# Function to annotate DE results
annotate_results <- function(tt) {
  genes <- tt$gene
  
  # Try to map gene symbols to Entrez IDs and descriptions
  tryCatch({
    # Map to Entrez ID
    entrez <- mapIds(org.Hs.eg.db, 
                     keys = genes, 
                     column = "ENTREZID", 
                     keytype = "SYMBOL",
                     multiVals = "first")
    
    # Map to gene name
    genename <- mapIds(org.Hs.eg.db, 
                       keys = genes, 
                       column = "GENENAME", 
                       keytype = "SYMBOL",
                       multiVals = "first")
    
    tt$entrez_id <- entrez[match(genes, names(entrez))]
    tt$gene_name <- genename[match(genes, names(genename))]
    
    message("✓ Gene annotations added")
  }, error = function(e) {
    message("! Could not add gene annotations: ", e$message)
    tt$entrez_id <- NA
    tt$gene_name <- NA
  })
  
  return(tt)
}

# Annotate each result set
for (contrast in names(de_results)) {
  de_results[[contrast]] <- annotate_results(de_results[[contrast]])
}

# =============================================================================
# 9. SAVE DE RESULTS
# =============================================================================

message("\n--- Saving DE Results ---\n")

# Save complete results for each contrast
for (contrast in names(de_results)) {
  # Full results
  file_full <- file.path(PATHS$deg_tables, paste0("DE_", contrast, "_full.csv"))
  write.csv(de_results[[contrast]], file_full, row.names = FALSE)
  message("✓ Saved: ", file_full)
  
  # Significant genes only
  sig_genes <- de_results[[contrast]][de_results[[contrast]]$significant, ]
  if (nrow(sig_genes) > 0) {
    file_sig <- file.path(PATHS$deg_tables, paste0("DE_", contrast, "_significant.csv"))
    write.csv(sig_genes, file_sig, row.names = FALSE)
    message("✓ Saved: ", file_sig)
  }
}

# Save summary table
summary_file <- file.path(PATHS$deg_tables, "DE_summary.csv")
write.csv(de_summary, summary_file, row.names = FALSE)
message("✓ Summary saved: ", summary_file)

# Save limma fit object for downstream analysis
fit_file <- file.path(PATHS$processed, "04_limma_fit.rds")
saveRDS(fit2, fit_file)
message("✓ limma fit object saved: ", fit_file)

# =============================================================================
# 10. DIAGNOSTIC PLOTS
# =============================================================================

message("\n--- Creating Diagnostic Plots ---\n")

# 10.1 Mean-Variance relationship (before and after eBayes)
pdf(file.path(PATHS$figures_de, "mean_variance_trend.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plotSA(fit2, main = "Mean-Variance Trend (after eBayes)")
dev.off()
message("✓ Mean-variance plot saved")

# 10.2 P-value distributions
pdf(file.path(PATHS$figures_de, "pvalue_distributions.pdf"), width = 12, height = 8)
par(mfrow = c(2, 2))
for (contrast in names(de_results)) {
  hist(de_results[[contrast]]$P.Value, 
       breaks = 50, 
       main = paste("P-value distribution:", contrast),
       xlab = "P-value",
       col = "lightblue",
       border = "white")
  # Add expected uniform distribution line
  abline(h = nrow(de_results[[contrast]]) / 50, col = "red", lty = 2)
}
dev.off()
message("✓ P-value distribution plots saved")

# 10.3 Volcano plots for each contrast
for (contrast in names(de_results)) {
  tt <- de_results[[contrast]]
  
  # Prepare data for volcano plot
  tt$neg_log10_pval <- -log10(tt$P.Value)
  tt$label <- ifelse(tt$significant, tt$gene, NA)
  
  # Only label top genes to avoid overcrowding
  top_genes <- head(tt$gene[tt$significant], 20)
  tt$label <- ifelse(tt$gene %in% top_genes, tt$gene, NA)
  
  # Create volcano plot
  p_volcano <- ggplot(tt, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), 
               linetype = "dashed", color = "gray40") +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 30,
                    segment.size = 0.2, segment.alpha = 0.5) +
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red")) +
    labs(
      title = paste("Volcano Plot:", contrast),
      subtitle = paste(sum(tt$significant), "significant genes"),
      x = "Log2 Fold Change",
      y = "-Log10 P-value",
      color = "Significant"
    ) +
    theme_publication() +
    theme(legend.position = "none")
  
  # Save
  ggsave(file.path(PATHS$figures_de, paste0("volcano_", contrast, ".pdf")),
         p_volcano, width = 10, height = 8, dpi = ANALYSIS_PARAMS$plot_dpi)
  ggsave(file.path(PATHS$figures_de, paste0("volcano_", contrast, ".png")),
         p_volcano, width = 10, height = 8, dpi = ANALYSIS_PARAMS$plot_dpi)
}
message("✓ Volcano plots saved")

# 10.4 MA plots
pdf(file.path(PATHS$figures_de, "MA_plots.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2))
for (i in seq_along(contrast_names)) {
  limma::plotMA(fit2, coef = i, status = decideTests(fit2)[, i],
                main = contrast_names[i],
                legend = FALSE)
}
dev.off()
message("✓ MA plots saved")

# =============================================================================
# 11. ANALYZE INTERACTION CONTRAST IN DETAIL
# =============================================================================

message("\n--- Interaction Contrast Analysis ---\n")

# The interaction contrast is scientifically most important
# It identifies TNBC-specific TP4 responses

interaction_results <- de_results[["Interaction"]]

message("Interaction contrast (TNBC-specific TP4 response):")
message("  Total significant: ", sum(interaction_results$significant))
message("  Up in TNBC: ", sum(interaction_results$significant & interaction_results$direction == "Up"))
message("  Down in TNBC: ", sum(interaction_results$significant & interaction_results$direction == "Down"))

# Interpretation guide
message("\nInterpretation of significant interaction genes:")
message("  Positive logFC: Gene is MORE upregulated by TP4 in TNBC vs HDF")
message("                  OR less downregulated in TNBC")
message("  Negative logFC: Gene is MORE downregulated by TP4 in TNBC vs HDF")
message("                  OR less upregulated in TNBC")

# Show top TNBC-specific upregulated genes
message("\nTop genes MORE upregulated by TP4 in TNBC (potential stress/death responses):")
up_interaction <- interaction_results %>%
  filter(significant, direction == "Up") %>%
  arrange(desc(logFC)) %>%
  head(20)
print(up_interaction[, c("gene", "logFC", "adj.P.Val", "gene_name")])

# Show top TNBC-specific downregulated genes  
message("\nTop genes MORE downregulated by TP4 in TNBC:")
down_interaction <- interaction_results %>%
  filter(significant, direction == "Down") %>%
  arrange(logFC) %>%
  head(20)
print(down_interaction[, c("gene", "logFC", "adj.P.Val", "gene_name")])

# =============================================================================
# 12. COMPARE CONTRASTS (VENN DIAGRAM DATA)
# =============================================================================

message("\n--- Comparing Contrasts ---\n")

# Get significant genes from each contrast
sig_lists <- lapply(de_results, function(x) {
  x$gene[x$significant]
})

# Calculate overlaps
overlap_matrix <- matrix(0, 
                         nrow = length(sig_lists), 
                         ncol = length(sig_lists),
                         dimnames = list(names(sig_lists), names(sig_lists)))

for (i in names(sig_lists)) {
  for (j in names(sig_lists)) {
    overlap_matrix[i, j] <- length(intersect(sig_lists[[i]], sig_lists[[j]]))
  }
}

message("Overlap matrix (significant genes):")
print(overlap_matrix)

# Genes unique to interaction
mda_genes <- sig_lists[["TP4_in_MDA"]]
hdf_genes <- sig_lists[["TP4_in_HDF"]]
interaction_genes <- sig_lists[["Interaction"]]

unique_to_interaction <- setdiff(interaction_genes, union(mda_genes, hdf_genes))
message("\nGenes unique to interaction contrast: ", length(unique_to_interaction))

# Common response genes (in both MDA and HDF)
common_response <- intersect(mda_genes, hdf_genes)
message("Common TP4 response genes (in both cell lines): ", length(common_response))

# MDA-specific response
mda_specific <- setdiff(mda_genes, hdf_genes)
message("MDA-MB-231 specific TP4 response genes: ", length(mda_specific))

# HDF-specific response
hdf_specific <- setdiff(hdf_genes, mda_genes)
message("HDF-specific TP4 response genes: ", length(hdf_specific))

# =============================================================================
# 13. CREATE DE OUTPUT OBJECT
# =============================================================================

message("\n--- Creating DE Output ---\n")

de_output <- list(
  # Results
  results = de_results,
  summary = de_summary,
  
  # limma fit
  fit = fit2,
  design = design,
  contrasts = contrast_matrix,
  
  # Gene lists
  gene_lists = list(
    TP4_in_MDA = sig_lists[["TP4_in_MDA"]],
    TP4_in_HDF = sig_lists[["TP4_in_HDF"]],
    Interaction = sig_lists[["Interaction"]],
    MDA_specific = mda_specific,
    HDF_specific = hdf_specific,
    Common = common_response,
    Unique_to_interaction = unique_to_interaction
  ),
  
  # Parameters
  parameters = list(
    fdr_threshold = fdr_threshold,
    log2fc_threshold = log2fc_threshold,
    analysis_level = analysis_level
  )
)

# Save complete DE output
output_file <- file.path(PATHS$processed, "04_de_output.rds")
saveRDS(de_output, output_file)
message("✓ Complete DE output saved: ", output_file)

# =============================================================================
# 14. GENERATE DE REPORT
# =============================================================================

de_report <- paste0(
  "================================================================================\n",
  "DIFFERENTIAL EXPRESSION ANALYSIS REPORT\n",
  "================================================================================\n\n",
  "Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "Dataset: GSE74764\n\n",
  
  "STATISTICAL MODEL\n",
  "-----------------\n",
  "Model: Cell-means model (Expression ~ 0 + Group)\n",
  "Groups: HDF_Mock, HDF_TP4, MDA-MB-231_Mock, MDA-MB-231_TP4\n",
  "Method: limma with empirical Bayes moderation\n",
  "Multiple testing: Benjamini-Hochberg FDR\n\n",
  
  "THRESHOLDS\n",
  "----------\n",
  "FDR: < ", fdr_threshold, "\n",
  "Log2 Fold-Change: > ", log2fc_threshold, " (", round(2^log2fc_threshold, 2), "-fold)\n\n",
  
  "RESULTS SUMMARY\n",
  "---------------\n",
  paste(capture.output(print(de_summary)), collapse = "\n"), "\n\n",
  
  "KEY FINDINGS\n",
  "------------\n",
  "1. TP4 effect in TNBC (MDA-MB-231): ", de_summary$Sig_FDR_FC[de_summary$Contrast == "TP4_in_MDA"], " DEGs\n",
  "2. TP4 effect in normal (HDF): ", de_summary$Sig_FDR_FC[de_summary$Contrast == "TP4_in_HDF"], " DEGs\n",
  "3. TNBC-specific TP4 response (Interaction): ", de_summary$Sig_FDR_FC[de_summary$Contrast == "Interaction"], " DEGs\n\n",
  
  "INTERPRETATION\n",
  "--------------\n",
  "The INTERACTION contrast identifies genes that respond differently to TP4\n",
  "in cancer cells vs normal fibroblasts. These genes are candidates for\n",
  "explaining TP4's selective anticancer activity.\n\n",
  
  "- Positive interaction logFC: Gene MORE induced by TP4 in TNBC\n",
  "- Negative interaction logFC: Gene MORE suppressed by TP4 in TNBC\n\n",
  
  "GENE OVERLAPS\n",
  "-------------\n",
  "MDA-specific response: ", length(mda_specific), " genes\n",
  "HDF-specific response: ", length(hdf_specific), " genes\n",
  "Common response: ", length(common_response), " genes\n\n",
  
  "CAVEATS\n",
  "-------\n",
  "1. Small sample size (n=2 per group) limits statistical power\n",
  "2. Results are hypothesis-generating, not definitive\n",
  "3. Validation in independent datasets is recommended\n\n",
  
  "OUTPUT FILES\n",
  "------------\n",
  "1. DE_[contrast]_full.csv - Complete results for each contrast\n",
  "2. DE_[contrast]_significant.csv - Significant genes only\n",
  "3. volcano_[contrast].pdf/png - Volcano plots\n",
  "4. MA_plots.pdf - MA plots\n",
  "5. pvalue_distributions.pdf - P-value histograms\n\n",
  
  "================================================================================\n"
)

# Save report
report_file <- file.path(PATHS$deg_tables, "DE_report.txt")
writeLines(de_report, report_file)
message("✓ DE report saved: ", report_file)

# Print report
cat(de_report)

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Differential expression complete. Proceed to 05_functional_analysis.R   ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

################################################################################
# END OF SCRIPT
################################################################################
