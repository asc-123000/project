################################################################################
# 01_data_download.R - GEO Data Acquisition
################################################################################
#
# Project: TP4-Induced Transcriptomic Changes in Triple-Negative Breast Cancer
# Dataset: GSE74764 (Agilent SurePrint G3 Human GE v2 8x60K)
# 
# Purpose:
#   Download the GSE74764 dataset from GEO, including:
#   - Series matrix file (processed expression data)
#   - Platform annotation (GPL16699)
#   - Sample metadata
#
# Scientific Rationale for Using Processed Data:
# ─────────────────────────────────────────────────
#   1. SMALL SAMPLE SIZE (n=2 per group):
#      With only 2 biological replicates per condition, any additional
#      technical variability introduced during raw data preprocessing
#      could disproportionately affect results. The GEO-processed data
#      has been normalized by the original authors who had full knowledge
#      of the experimental design.
#
#   2. PLATFORM-SPECIFIC EXPERTISE:
#      Agilent two-color arrays require specialized background correction
#      and normalization. The original authors applied appropriate methods
#      for their specific experimental protocol.
#
#   3. REPRODUCIBILITY:
#      Using the same processed data as the original publication allows
#      for more direct comparison with their findings.
#
#   4. QC VERIFICATION:
#      We will perform thorough QC to verify the normalization quality
#      and identify any issues with the processed data.
#
# Note:
#   If QC reveals problems with the processed data, we can revisit
#   this decision and implement raw data preprocessing.
#
# Author: Computational Biology Pipeline
# Date: 2026-01-19
#
################################################################################

# =============================================================================
# 1. SETUP
# =============================================================================

# Clear workspace (optional - comment out if running as part of pipeline)
# rm(list = ls())

# Source the setup script to load packages and paths
# Adjust path as needed for your system
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
if (length(script_dir) == 0 || script_dir == "") {
  script_dir <- "scripts"  # Fallback for non-RStudio execution
}
source(file.path(dirname(script_dir), "scripts", "00_setup.R"))

# Alternatively, if running from project root:
# source("scripts/00_setup.R")

# Load required packages
library(GEOquery)
library(Biobase)
library(tidyverse)

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Script 01: GEO Data Download                                            ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

# =============================================================================
# 2. DOWNLOAD GSE74764 FROM GEO
# =============================================================================

message("Downloading GSE74764 from GEO...")
message("This may take a few minutes depending on internet connection.\n")

# Download the GEO dataset
# getGEO() returns a list of ExpressionSet objects
# GSE74764 should have one element corresponding to the GPL16699 platform

gse <- tryCatch({
  GEOquery::getGEO(
    GEO = ANALYSIS_PARAMS$geo_accession,
    destdir = PATHS$raw,
    GSEMatrix = TRUE,      # Get the processed expression matrix
    AnnotGPL = TRUE,       # Get annotation from GEO
    getGPL = TRUE          # Download platform annotation
  )
}, error = function(e) {
  message("Error downloading from GEO: ", e$message)
  message("\nTrying alternative download method...")
  
  # Alternative: download directly
  GEOquery::getGEO(
    GEO = ANALYSIS_PARAMS$geo_accession,
    destdir = PATHS$raw,
    GSEMatrix = TRUE,
    AnnotGPL = FALSE,
    getGPL = TRUE
  )
})

# The result is a list - extract the ExpressionSet
if (is.list(gse) && length(gse) >= 1) {
  eset <- gse[[1]]
  message("✓ Successfully downloaded GSE74764")
  message("  Platform: ", annotation(eset))
  message("  Samples: ", ncol(eset))
  message("  Features: ", nrow(eset))
} else {
  stop("Failed to retrieve ExpressionSet from GEO")
}

# =============================================================================
# 3. EXPLORE AND VERIFY DATASET STRUCTURE
# =============================================================================

message("\n--- Dataset Overview ---\n")

# Basic information
message("ExpressionSet dimensions:")
message("  Probes (features): ", nrow(eset))
message("  Samples: ", ncol(eset))

# Check expression data
expr_matrix <- Biobase::exprs(eset)
message("\nExpression matrix summary:")
message("  Range: [", round(min(expr_matrix, na.rm = TRUE), 2), ", ", 
        round(max(expr_matrix, na.rm = TRUE), 2), "]")
message("  Median: ", round(median(expr_matrix, na.rm = TRUE), 2))
message("  Missing values: ", sum(is.na(expr_matrix)))

# Check if data appears to be log-transformed
# Log2 data typically has values in range 0-16 for microarray
if (max(expr_matrix, na.rm = TRUE) < 20 && min(expr_matrix, na.rm = TRUE) > -5) {
  message("  Data appears to be log2-transformed ✓")
  is_log2 <- TRUE
} else {
  message("  Data may NOT be log2-transformed - will need transformation")
  is_log2 <- FALSE
}

# =============================================================================
# 4. EXTRACT AND PARSE SAMPLE METADATA
# =============================================================================

message("\n--- Sample Metadata ---\n")

# Get phenotype data (sample metadata)
pheno_data <- Biobase::pData(eset)

# Display available columns
message("Available metadata columns:")
print(colnames(pheno_data))

# Extract key information from sample characteristics
# The structure varies by dataset, so we need to explore

# Display sample titles and characteristics
message("\nSample information:")
for (i in 1:ncol(eset)) {
  message("  ", i, ": ", pheno_data$title[i])
}

# =============================================================================
# 5. CREATE CLEAN SAMPLE METADATA
# =============================================================================

message("\n--- Creating Clean Sample Metadata ---\n")

# Parse sample information from titles or characteristics
# Expected format: cell line + treatment information

# Create a clean metadata data frame
sample_metadata <- data.frame(
  sample_id = colnames(eset),
  geo_accession = pheno_data$geo_accession,
  title = pheno_data$title,
  stringsAsFactors = FALSE
)

# Parse cell line and treatment from titles
# Typical format: "MDA-MB-231_TP4_rep1" or similar
# We need to examine the actual titles to parse correctly

# Function to extract cell line from title
extract_cell_line <- function(title) {
  if (grepl("MDA|231", title, ignore.case = TRUE)) {
    return("MDA-MB-231")
  } else if (grepl("HDF|fibroblast", title, ignore.case = TRUE)) {
    return("HDF")
  } else {
    return(NA)
  }
}

# Function to extract treatment from title
extract_treatment <- function(title) {
  if (grepl("TP4|treated|piscidin", title, ignore.case = TRUE) && 
      !grepl("mock|control|untreated", title, ignore.case = TRUE)) {
    return("TP4")
  } else if (grepl("mock|control|untreated", title, ignore.case = TRUE)) {
    return("Mock")
  } else {
    # Check characteristics columns if available
    return(NA)
  }
}

# Apply parsing functions
sample_metadata$cell_line <- sapply(sample_metadata$title, extract_cell_line)
sample_metadata$treatment <- sapply(sample_metadata$title, extract_treatment)

# If parsing from titles failed, try characteristics columns
if (any(is.na(sample_metadata$cell_line)) || any(is.na(sample_metadata$treatment))) {
  message("Parsing from titles incomplete. Checking characteristics columns...")
  
  # Look for characteristics columns
  char_cols <- grep("characteristics", colnames(pheno_data), value = TRUE, ignore.case = TRUE)
  
  if (length(char_cols) > 0) {
    for (col in char_cols) {
      message("  Checking: ", col)
      message("    Values: ", paste(unique(pheno_data[[col]]), collapse = "; "))
    }
    
    # Re-parse using characteristics
    for (col in char_cols) {
      values <- pheno_data[[col]]
      
      # Check for cell line info
      if (any(grepl("cell", col, ignore.case = TRUE)) || 
          any(grepl("MDA|HDF|231|fibroblast", values, ignore.case = TRUE))) {
        sample_metadata$cell_line <- sapply(values, function(v) {
          if (grepl("MDA|231", v, ignore.case = TRUE)) "MDA-MB-231"
          else if (grepl("HDF|fibroblast", v, ignore.case = TRUE)) "HDF"
          else NA
        })
      }
      
      # Check for treatment info
      if (any(grepl("treatment|agent", col, ignore.case = TRUE)) ||
          any(grepl("TP4|mock|control", values, ignore.case = TRUE))) {
        sample_metadata$treatment <- sapply(values, function(v) {
          if (grepl("TP4|piscidin", v, ignore.case = TRUE) && 
              !grepl("mock|control", v, ignore.case = TRUE)) "TP4"
          else if (grepl("mock|control|untreated", v, ignore.case = TRUE)) "Mock"
          else NA
        })
      }
    }
  }
}

# Create group factor (combination of cell line and treatment)
sample_metadata$group <- paste(sample_metadata$cell_line, sample_metadata$treatment, sep = "_")

# Add replicate information
sample_metadata <- sample_metadata %>%
  group_by(group) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  as.data.frame()

# Verify we have the expected 2x2 factorial design with 2 replicates each
message("\nSample metadata summary:")
print(table(sample_metadata$cell_line, sample_metadata$treatment))

# Verify expected sample count
expected_samples <- 8
if (nrow(sample_metadata) != expected_samples) {
  warning("Expected ", expected_samples, " samples but found ", nrow(sample_metadata))
}

# Check for missing values
if (any(is.na(sample_metadata$cell_line)) || any(is.na(sample_metadata$treatment))) {
  warning("Some samples could not be assigned to cell line or treatment groups!")
  message("Problematic samples:")
  print(sample_metadata[is.na(sample_metadata$cell_line) | is.na(sample_metadata$treatment), ])
  
  # Manual assignment may be needed - display all info for user review
  message("\nFull phenotype data for manual review:")
  print(pheno_data[, 1:min(10, ncol(pheno_data))])
}

# Display final metadata
message("\nFinal sample metadata:")
print(sample_metadata[, c("sample_id", "cell_line", "treatment", "group", "replicate")])

# =============================================================================
# 6. EXTRACT PROBE ANNOTATION
# =============================================================================

message("\n--- Probe Annotation ---\n")

# Get feature data (probe annotations)
feature_data <- Biobase::fData(eset)

message("Available annotation columns:")
print(colnames(feature_data))

# Identify key annotation columns
# Common column names for gene symbols and IDs
gene_symbol_cols <- c("Gene Symbol", "Gene.Symbol", "GENE_SYMBOL", "gene_symbol",
                      "Symbol", "GENE", "gene", "Gene symbol", "SYMBOL")
gene_id_cols <- c("Gene ID", "ENTREZ_GENE_ID", "Entrez_Gene_ID", "GENE_ID",
                  "gene_id", "EntrezID", "ENTREZID")
gene_name_cols <- c("Gene Name", "Gene_Name", "GENE_NAME", "gene_name",
                    "Gene title", "GENE_TITLE", "Description")

# Find matching columns
symbol_col <- intersect(gene_symbol_cols, colnames(feature_data))[1]
id_col <- intersect(gene_id_cols, colnames(feature_data))[1]
name_col <- intersect(gene_name_cols, colnames(feature_data))[1]

message("\nIdentified annotation columns:")
message("  Gene Symbol: ", ifelse(is.na(symbol_col), "NOT FOUND", symbol_col))
message("  Gene ID: ", ifelse(is.na(id_col), "NOT FOUND", id_col))
message("  Gene Name: ", ifelse(is.na(name_col), "NOT FOUND", name_col))

# Create clean annotation data frame
probe_annotation <- data.frame(
  probe_id = rownames(feature_data),
  stringsAsFactors = FALSE
)

if (!is.na(symbol_col)) {
  probe_annotation$gene_symbol <- feature_data[[symbol_col]]
}
if (!is.na(id_col)) {
  probe_annotation$entrez_id <- feature_data[[id_col]]
}
if (!is.na(name_col)) {
  probe_annotation$gene_name <- feature_data[[name_col]]
}

# Statistics on annotation coverage
if ("gene_symbol" %in% colnames(probe_annotation)) {
  annotated_probes <- sum(!is.na(probe_annotation$gene_symbol) & 
                          probe_annotation$gene_symbol != "" &
                          probe_annotation$gene_symbol != "---")
  message("\nAnnotation coverage:")
  message("  Total probes: ", nrow(probe_annotation))
  message("  Annotated probes: ", annotated_probes, 
          " (", round(100 * annotated_probes / nrow(probe_annotation), 1), "%)")
}

# =============================================================================
# 7. SAVE DOWNLOADED DATA
# =============================================================================

message("\n--- Saving Data ---\n")

# Save the ExpressionSet object
eset_file <- file.path(PATHS$raw, "GSE74764_eset.rds")
saveRDS(eset, eset_file)
message("✓ ExpressionSet saved: ", eset_file)

# Save expression matrix as CSV (for external tools)
expr_file <- file.path(PATHS$raw, "GSE74764_expression_matrix.csv")
expr_df <- as.data.frame(expr_matrix)
expr_df$probe_id <- rownames(expr_df)
expr_df <- expr_df[, c("probe_id", setdiff(colnames(expr_df), "probe_id"))]
write.csv(expr_df, expr_file, row.names = FALSE)
message("✓ Expression matrix saved: ", expr_file)

# Save sample metadata
metadata_file <- file.path(PATHS$raw, "sample_metadata.csv")
write.csv(sample_metadata, metadata_file, row.names = FALSE)
message("✓ Sample metadata saved: ", metadata_file)

# Save probe annotation
annotation_file <- file.path(PATHS$raw, "probe_annotation.csv")
write.csv(probe_annotation, annotation_file, row.names = FALSE)
message("✓ Probe annotation saved: ", annotation_file)

# Save full phenotype data for reference
pheno_file <- file.path(PATHS$raw, "GSE74764_phenotype_data.csv")
write.csv(pheno_data, pheno_file, row.names = TRUE)
message("✓ Full phenotype data saved: ", pheno_file)

# =============================================================================
# 8. DATA INTEGRITY CHECKS
# =============================================================================

message("\n--- Data Integrity Checks ---\n")

# Check 1: Expression matrix dimensions match metadata
check1 <- ncol(expr_matrix) == nrow(sample_metadata)
message("Check 1 - Sample counts match: ", ifelse(check1, "PASS ✓", "FAIL ✗"))

# Check 2: Sample IDs match between expression matrix and metadata
check2 <- all(colnames(expr_matrix) == sample_metadata$sample_id)
message("Check 2 - Sample IDs match: ", ifelse(check2, "PASS ✓", "FAIL ✗"))

# Check 3: All groups have replicates
group_counts <- table(sample_metadata$group)
check3 <- all(group_counts >= 2)
message("Check 3 - All groups have ≥2 replicates: ", ifelse(check3, "PASS ✓", "FAIL ✗"))

# Check 4: No missing expression values (or acceptable level)
missing_pct <- 100 * sum(is.na(expr_matrix)) / length(expr_matrix)
check4 <- missing_pct < 5  # Less than 5% missing
message("Check 4 - Missing values < 5%: ", ifelse(check4, "PASS ✓", "FAIL ✗"),
        " (", round(missing_pct, 2), "%)")

# Overall status
all_checks_pass <- all(check1, check2, check3, check4)

if (all_checks_pass) {
  message("\n✓ All data integrity checks passed!")
} else {
  warning("\n✗ Some data integrity checks failed - please review")
}

# =============================================================================
# 9. CREATE DATA SUMMARY REPORT
# =============================================================================

message("\n--- Creating Data Summary Report ---\n")

summary_report <- paste0(
  "================================================================================\n",
  "GSE74764 DATA DOWNLOAD SUMMARY\n",
  "================================================================================\n\n",
  "Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
  "DATASET INFORMATION\n",
  "-------------------\n",
  "GEO Accession: ", ANALYSIS_PARAMS$geo_accession, "\n",
  "Platform: ", annotation(eset), "\n",
  "Organism: Homo sapiens\n\n",
  "SAMPLE SUMMARY\n",
  "--------------\n",
  "Total samples: ", ncol(eset), "\n",
  "Cell lines: ", paste(unique(sample_metadata$cell_line), collapse = ", "), "\n",
  "Treatments: ", paste(unique(sample_metadata$treatment), collapse = ", "), "\n\n",
  "EXPRESSION DATA\n",
  "---------------\n",
  "Total probes: ", nrow(eset), "\n",
  "Data range: [", round(min(expr_matrix, na.rm = TRUE), 2), ", ", 
  round(max(expr_matrix, na.rm = TRUE), 2), "]\n",
  "Log2-transformed: ", ifelse(is_log2, "Yes", "No"), "\n",
  "Missing values: ", sum(is.na(expr_matrix)), " (", round(missing_pct, 2), "%)\n\n",
  "ANNOTATION COVERAGE\n",
  "-------------------\n",
  "Probes with gene symbols: ", 
  ifelse("gene_symbol" %in% colnames(probe_annotation),
         paste0(annotated_probes, " (", round(100 * annotated_probes / nrow(probe_annotation), 1), "%)"),
         "N/A"), "\n\n",
  "OUTPUT FILES\n",
  "------------\n",
  "1. ", eset_file, "\n",
  "2. ", expr_file, "\n",
  "3. ", metadata_file, "\n",
  "4. ", annotation_file, "\n",
  "5. ", pheno_file, "\n\n",
  "DATA INTEGRITY\n",
  "--------------\n",
  "All checks passed: ", ifelse(all_checks_pass, "Yes", "No"), "\n\n",
  "================================================================================\n"
)

# Save summary report
report_file <- file.path(PATHS$raw, "download_summary.txt")
writeLines(summary_report, report_file)
message("✓ Summary report saved: ", report_file)

# Print summary
cat(summary_report)

# =============================================================================
# 10. PREPARE OBJECTS FOR NEXT SCRIPT
# =============================================================================

# Create a list of objects to pass to the next script
download_output <- list(
  eset = eset,
  expression_matrix = expr_matrix,
  sample_metadata = sample_metadata,
  probe_annotation = probe_annotation,
  is_log2 = is_log2,
  files = list(
    eset = eset_file,
    expression = expr_file,
    metadata = metadata_file,
    annotation = annotation_file
  )
)

# Save for next script
output_file <- file.path(PATHS$processed, "01_download_output.rds")
saveRDS(download_output, output_file)
message("\n✓ Download output saved for next script: ", output_file)

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Data download complete. Proceed to 02_preprocessing.R                   ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

################################################################################
# END OF SCRIPT
################################################################################
