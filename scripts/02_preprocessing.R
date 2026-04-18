################################################################################
# 02_preprocessing.R - Data Preprocessing and Normalization Verification
################################################################################
#
# Project: TP4-Induced Transcriptomic Changes in Triple-Negative Breast Cancer
# Dataset: GSE74764 (Agilent SurePrint G3 Human GE v2 8x60K)
# 
# Purpose:
#   Process the GEO series matrix data:
#   1. Verify normalization quality of pre-processed data
#   2. Handle missing values (if any)
#   3. Filter low-expression probes
#   4. Map probes to genes and handle multi-mapping
#   5. Prepare analysis-ready expression matrix
#
# Scientific Rationale:
# ─────────────────────────────────────────────────
#   We use the GEO-processed data because:
#   
#   1. The small sample size (n=2 per group) makes our analysis sensitive
#      to any additional technical variability. Raw data preprocessing
#      would require choices (background correction method, normalization
#      algorithm) that could introduce variance.
#   
#   2. The original authors applied Agilent-specific processing appropriate
#      for their two-color array protocol.
#   
#   3. This script verifies the quality of the pre-processed data and
#      prepares it for differential expression analysis.
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
library(Biobase)
library(tidyverse)

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Script 02: Data Preprocessing                                           ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

# =============================================================================
# 2. LOAD DOWNLOADED DATA
# =============================================================================

message("--- Loading Downloaded Data ---\n")

# Load the output from the previous script
download_output <- readRDS(file.path(PATHS$processed, "01_download_output.rds"))

eset <- download_output$eset
expr_matrix <- download_output$expression_matrix
sample_metadata <- download_output$sample_metadata
probe_annotation <- download_output$probe_annotation
is_log2 <- download_output$is_log2

message("✓ Loaded expression data: ", nrow(expr_matrix), " probes × ", ncol(expr_matrix), " samples")
message("✓ Sample groups: ", paste(unique(sample_metadata$group), collapse = ", "))

# =============================================================================
# 3. LOG2 TRANSFORMATION (IF NEEDED)
# =============================================================================

message("\n--- Checking Log2 Transformation ---\n")

# Check data distribution to confirm log2 status
data_range <- range(expr_matrix, na.rm = TRUE)
data_median <- median(expr_matrix, na.rm = TRUE)

message("Expression data range: [", round(data_range[1], 2), ", ", round(data_range[2], 2), "]")
message("Expression data median: ", round(data_median, 2))

# Log2 microarray data typically has:
# - Range approximately 2-16
# - Median around 6-10
if (!is_log2 && data_range[2] > 100) {
  message("\nData appears to be on linear scale. Applying log2 transformation...")
  
  # Add small offset to avoid log(0)
  min_val <- min(expr_matrix[expr_matrix > 0], na.rm = TRUE)
  offset <- min_val / 2
  
  expr_matrix <- log2(expr_matrix + offset)
  
  message("✓ Log2 transformation applied (offset = ", round(offset, 4), ")")
  message("  New range: [", round(min(expr_matrix, na.rm = TRUE), 2), ", ", 
          round(max(expr_matrix, na.rm = TRUE), 2), "]")
} else {
  message("✓ Data is already log2-transformed")
}

# =============================================================================
# 4. HANDLE MISSING VALUES
# =============================================================================

message("\n--- Handling Missing Values ---\n")

# Count missing values
n_missing <- sum(is.na(expr_matrix))
pct_missing <- 100 * n_missing / length(expr_matrix)

message("Missing values: ", n_missing, " (", round(pct_missing, 2), "%)")

if (n_missing > 0) {
  # Identify probes with missing values
  probes_with_na <- rowSums(is.na(expr_matrix))
  probes_any_na <- sum(probes_with_na > 0)
  probes_all_na <- sum(probes_with_na == ncol(expr_matrix))
  
  message("  Probes with any missing: ", probes_any_na)
  message("  Probes with all missing: ", probes_all_na)
  
  # Remove probes that are entirely missing
  if (probes_all_na > 0) {
    expr_matrix <- expr_matrix[probes_with_na < ncol(expr_matrix), ]
    message("  Removed ", probes_all_na, " probes with all missing values")
  }
  
  # For probes with partial missing, impute with row median
  # This is a conservative approach for microarray data
  remaining_na <- sum(is.na(expr_matrix))
  if (remaining_na > 0) {
    message("  Imputing ", remaining_na, " remaining missing values with row medians...")
    
    for (i in which(rowSums(is.na(expr_matrix)) > 0)) {
      row_median <- median(expr_matrix[i, ], na.rm = TRUE)
      expr_matrix[i, is.na(expr_matrix[i, ])] <- row_median
    }
    
    message("  ✓ Imputation complete")
  }
} else {
  message("✓ No missing values to handle")
}

# =============================================================================
# 5. VERIFY NORMALIZATION QUALITY
# =============================================================================

message("\n--- Verifying Normalization Quality ---\n")

# Calculate per-sample statistics
sample_stats <- data.frame(
  sample_id = colnames(expr_matrix),
  mean = colMeans(expr_matrix),
  median = apply(expr_matrix, 2, median),
  sd = apply(expr_matrix, 2, sd),
  iqr = apply(expr_matrix, 2, IQR),
  min = apply(expr_matrix, 2, min),
  max = apply(expr_matrix, 2, max)
)

# Add metadata
sample_stats <- merge(sample_stats, sample_metadata[, c("sample_id", "cell_line", "treatment", "group")],
                      by = "sample_id")

message("Per-sample statistics:")
print(sample_stats[, c("sample_id", "group", "mean", "median", "sd")])

# Check for normalization quality
# Well-normalized data should have similar medians and IQRs
median_cv <- sd(sample_stats$median) / mean(sample_stats$median) * 100
iqr_cv <- sd(sample_stats$iqr) / mean(sample_stats$iqr) * 100

message("\nNormalization quality metrics:")
message("  Median CV across samples: ", round(median_cv, 2), "%")
message("  IQR CV across samples: ", round(iqr_cv, 2), "%")

if (median_cv < 10 && iqr_cv < 15) {
  message("  ✓ Data appears well-normalized")
  normalization_ok <- TRUE
} else {
  message("  ! Normalization may be suboptimal - will proceed with caution")
  normalization_ok <- FALSE
}

# Save sample statistics
stats_file <- file.path(PATHS$processed, "sample_statistics.csv")
write.csv(sample_stats, stats_file, row.names = FALSE)
message("\n✓ Sample statistics saved: ", stats_file)

# =============================================================================
# 6. FILTER LOW-EXPRESSION PROBES
# =============================================================================

message("\n--- Filtering Low-Expression Probes ---\n")

# Rationale: Probes with very low expression across all samples contribute
# noise rather than signal. Filtering improves statistical power by reducing
# the multiple testing burden.

# Calculate detection metrics
# A probe should be "expressed" in at least some samples
# We define "expressed" as being above a threshold (e.g., median of all values)

initial_probes <- nrow(expr_matrix)
global_median <- median(expr_matrix)

message("Initial probes: ", initial_probes)
message("Global median expression: ", round(global_median, 2))

# Calculate the proportion of samples where each probe is above the threshold
# Threshold: 10th percentile of all expression values (lenient)
expression_threshold <- quantile(expr_matrix, 0.10)
message("Expression threshold (10th percentile): ", round(expression_threshold, 2))

# Count samples above threshold for each probe
samples_above <- rowSums(expr_matrix > expression_threshold)

# Require expression in at least 2 samples (minimum for any group)
# This is lenient because we have small sample sizes
min_samples <- 2
keep_probes <- samples_above >= min_samples

expr_filtered <- expr_matrix[keep_probes, ]

message("Probes above threshold in ≥", min_samples, " samples: ", sum(keep_probes))
message("Probes removed: ", sum(!keep_probes), 
        " (", round(100 * sum(!keep_probes) / initial_probes, 1), "%)")

# Also remove probes with very low variance (likely noise or control probes)
probe_variance <- apply(expr_filtered, 1, var)
variance_threshold <- quantile(probe_variance, 0.05)  # Remove bottom 5%
keep_variance <- probe_variance > variance_threshold

expr_filtered <- expr_filtered[keep_variance, ]

message("Probes with low variance removed: ", sum(!keep_variance))
message("Final probe count: ", nrow(expr_filtered))

# Update probe annotation
probe_annotation_filtered <- probe_annotation[probe_annotation$probe_id %in% rownames(expr_filtered), ]

# =============================================================================
# 7. PROBE-TO-GENE MAPPING
# =============================================================================

message("\n--- Probe-to-Gene Mapping ---\n")

# Check annotation coverage
if ("gene_symbol" %in% colnames(probe_annotation_filtered)) {
  # Clean gene symbols
  probe_annotation_filtered$gene_symbol_clean <- probe_annotation_filtered$gene_symbol
  
  # Remove empty, NA, or placeholder values
  invalid_symbols <- c("", "---", "-", "NA", "N/A", "null")
  probe_annotation_filtered$gene_symbol_clean[
    probe_annotation_filtered$gene_symbol_clean %in% invalid_symbols |
    is.na(probe_annotation_filtered$gene_symbol_clean)
  ] <- NA
  
  annotated_probes <- sum(!is.na(probe_annotation_filtered$gene_symbol_clean))
  message("Probes with valid gene symbols: ", annotated_probes, 
          " (", round(100 * annotated_probes / nrow(probe_annotation_filtered), 1), "%)")
} else {
  message("! Gene symbol column not found in annotation")
  message("  Will attempt to map using probe IDs...")
  
  # Add placeholder
  probe_annotation_filtered$gene_symbol_clean <- NA
}

# =============================================================================
# 8. HANDLE MULTI-MAPPING PROBES
# =============================================================================

message("\n--- Handling Multi-Mapping Probes ---\n")

# Some probes map to multiple genes (e.g., "GENE1 /// GENE2")
# Strategy: For probes mapping to multiple genes, we'll keep the first gene
# This is a common approach; alternatives include averaging across genes

if ("gene_symbol_clean" %in% colnames(probe_annotation_filtered)) {
  # Check for multi-mapping indicators
  multi_map_indicator <- " /// "
  
  multi_mapping <- grepl(multi_map_indicator, probe_annotation_filtered$gene_symbol_clean, fixed = TRUE)
  message("Probes mapping to multiple genes: ", sum(multi_mapping, na.rm = TRUE))
  
  if (sum(multi_mapping, na.rm = TRUE) > 0) {
    # Take the first gene symbol
    probe_annotation_filtered$gene_symbol_primary <- sapply(
      probe_annotation_filtered$gene_symbol_clean,
      function(x) {
        if (is.na(x)) return(NA)
        strsplit(x, " /// ", fixed = TRUE)[[1]][1]
      }
    )
    message("  → Using first gene symbol for multi-mapping probes")
  } else {
    probe_annotation_filtered$gene_symbol_primary <- probe_annotation_filtered$gene_symbol_clean
  }
  
  # Report unique genes
  unique_genes <- length(unique(na.omit(probe_annotation_filtered$gene_symbol_primary)))
  message("Unique gene symbols: ", unique_genes)
}

# =============================================================================
# 9. CREATE GENE-LEVEL EXPRESSION MATRIX
# =============================================================================

message("\n--- Creating Gene-Level Expression Matrix ---\n")

# When multiple probes map to the same gene, we need to collapse them
# Strategy: Use the probe with highest average expression
# This is biologically motivated: the highest-expressed probe likely
# represents the most relevant transcript

if ("gene_symbol_primary" %in% colnames(probe_annotation_filtered)) {
  
  # Add probe expression mean
  probe_annotation_filtered$mean_expr <- rowMeans(expr_filtered[probe_annotation_filtered$probe_id, ])
  
  # For annotated probes, select the best probe per gene
  annotated <- probe_annotation_filtered[!is.na(probe_annotation_filtered$gene_symbol_primary), ]
  
  best_probes <- annotated %>%
    group_by(gene_symbol_primary) %>%
    slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    pull(probe_id)
  
  message("Selecting best probe per gene (highest mean expression)")
  message("  Annotated probes: ", nrow(annotated))
  message("  Unique genes: ", length(best_probes))
  
  # Create gene-level expression matrix
  expr_gene <- expr_filtered[best_probes, ]
  
  # Replace rownames with gene symbols
  gene_mapping <- setNames(
    annotated$gene_symbol_primary[match(best_probes, annotated$probe_id)],
    best_probes
  )
  rownames(expr_gene) <- gene_mapping[rownames(expr_gene)]
  
  message("Gene-level expression matrix: ", nrow(expr_gene), " genes × ", ncol(expr_gene), " samples")
  
  # Also keep probe-level matrix for reference
  expr_probe <- expr_filtered
  
} else {
  message("! Cannot create gene-level matrix without gene annotations")
  message("  Using probe-level data for analysis")
  
  expr_gene <- NULL
  expr_probe <- expr_filtered
}

# =============================================================================
# 10. FINAL DATA PREPARATION
# =============================================================================

message("\n--- Final Data Preparation ---\n")

# Ensure sample metadata is in the same order as expression matrix
sample_metadata <- sample_metadata[match(colnames(expr_filtered), sample_metadata$sample_id), ]

# Verify order
stopifnot(all(colnames(expr_filtered) == sample_metadata$sample_id))
message("✓ Sample order verified")

# Create factor variables for analysis
sample_metadata$cell_line <- factor(sample_metadata$cell_line, 
                                     levels = c("HDF", "MDA-MB-231"))
sample_metadata$treatment <- factor(sample_metadata$treatment, 
                                     levels = c("Mock", "TP4"))
sample_metadata$group <- factor(sample_metadata$group,
                                 levels = c("HDF_Mock", "HDF_TP4", 
                                           "MDA-MB-231_Mock", "MDA-MB-231_TP4"))

message("Factor levels:")
message("  Cell line: ", paste(levels(sample_metadata$cell_line), collapse = " < "))
message("  Treatment: ", paste(levels(sample_metadata$treatment), collapse = " < "))
message("  Group: ", paste(levels(sample_metadata$group), collapse = ", "))

# =============================================================================
# 11. SAVE PROCESSED DATA
# =============================================================================

message("\n--- Saving Processed Data ---\n")

# Save probe-level expression matrix
probe_file <- file.path(PATHS$processed, "expression_probe_level.rds")
saveRDS(expr_probe, probe_file)
message("✓ Probe-level expression saved: ", probe_file)

# Save gene-level expression matrix
if (!is.null(expr_gene)) {
  gene_file <- file.path(PATHS$processed, "expression_gene_level.rds")
  saveRDS(expr_gene, gene_file)
  message("✓ Gene-level expression saved: ", gene_file)
  
  # Also save as CSV for external use
  gene_csv <- file.path(PATHS$processed, "expression_gene_level.csv")
  expr_gene_df <- as.data.frame(expr_gene)
  expr_gene_df$gene <- rownames(expr_gene_df)
  expr_gene_df <- expr_gene_df[, c("gene", setdiff(colnames(expr_gene_df), "gene"))]
  write.csv(expr_gene_df, gene_csv, row.names = FALSE)
  message("✓ Gene-level expression (CSV) saved: ", gene_csv)
}

# Save sample metadata
metadata_file <- file.path(PATHS$processed, "sample_metadata_processed.rds")
saveRDS(sample_metadata, metadata_file)
message("✓ Sample metadata saved: ", metadata_file)

# Save probe annotation
annotation_file <- file.path(PATHS$processed, "probe_annotation_filtered.rds")
saveRDS(probe_annotation_filtered, annotation_file)
message("✓ Probe annotation saved: ", annotation_file)

# =============================================================================
# 12. CREATE PREPROCESSING SUMMARY
# =============================================================================

message("\n--- Preprocessing Summary ---\n")

summary_list <- list(
  # Data matrices
  expression_probe = expr_probe,
  expression_gene = expr_gene,
  
  # Metadata
  sample_metadata = sample_metadata,
  probe_annotation = probe_annotation_filtered,
  
  # Processing parameters
  parameters = list(
    log2_transformed = TRUE,
    expression_threshold = expression_threshold,
    min_samples_expressed = min_samples,
    variance_threshold = variance_threshold,
    multi_mapping_strategy = "first_gene",
    probe_selection_strategy = "highest_mean_expression"
  ),
  
  # Summary statistics
  summary = list(
    initial_probes = initial_probes,
    filtered_probes = nrow(expr_probe),
    unique_genes = ifelse(!is.null(expr_gene), nrow(expr_gene), NA),
    samples = ncol(expr_probe),
    normalization_quality = normalization_ok
  ),
  
  # Sample statistics
  sample_stats = sample_stats
)

# Save complete preprocessing output
output_file <- file.path(PATHS$processed, "02_preprocessing_output.rds")
saveRDS(summary_list, output_file)
message("✓ Preprocessing output saved: ", output_file)

# Print summary
cat("\n")
cat("================================================================================\n")
cat("PREPROCESSING SUMMARY\n")
cat("================================================================================\n\n")
cat("INPUT:\n")
cat("  Initial probes: ", initial_probes, "\n")
cat("  Samples: ", ncol(expr_matrix), "\n\n")
cat("FILTERING:\n")
cat("  Low expression filter: ", sum(!keep_probes), " probes removed\n")
cat("  Low variance filter: ", sum(!keep_variance), " probes removed\n\n")
cat("OUTPUT:\n")
cat("  Probe-level features: ", nrow(expr_probe), "\n")
cat("  Gene-level features: ", ifelse(!is.null(expr_gene), nrow(expr_gene), "N/A"), "\n")
cat("  Samples: ", ncol(expr_probe), "\n\n")
cat("QUALITY:\n")
cat("  Normalization: ", ifelse(normalization_ok, "Good", "Suboptimal"), "\n")
cat("  Missing values: ", n_missing, " (imputed)\n")
cat("================================================================================\n")

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Preprocessing complete. Proceed to 03_quality_control.R                 ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

################################################################################
# END OF SCRIPT
################################################################################
