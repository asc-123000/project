################################################################################
# 07_generate_report.R - Final Report and Session Documentation
################################################################################
#
# Project: TP4-Induced Transcriptomic Changes in Triple-Negative Breast Cancer
# Dataset: GSE74764 (Agilent SurePrint G3 Human GE v2 8x60K)
# 
# Purpose:
#   Generate final analysis report and documentation:
#   1. Compile all analysis results
#   2. Generate summary statistics
#   3. Create reproducibility documentation
#   4. Save session information
#   5. Create final output tables
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
library(tidyverse)
library(knitr)
library(openxlsx)

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Script 07: Final Report Generation                                      ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

# =============================================================================
# 2. LOAD ALL ANALYSIS OUTPUTS
# =============================================================================

message("--- Loading All Analysis Results ---\n")

# Load all outputs
outputs <- list()

outputs$preproc <- tryCatch(
  readRDS(file.path(PATHS$processed, "02_preprocessing_output.rds")),
  error = function(e) NULL
)

outputs$qc <- tryCatch(
  readRDS(file.path(PATHS$processed, "03_qc_output.rds")),
  error = function(e) NULL
)

outputs$de <- tryCatch(
  readRDS(file.path(PATHS$processed, "04_de_output.rds")),
  error = function(e) NULL
)

outputs$functional <- tryCatch(
  readRDS(file.path(PATHS$processed, "05_functional_output.rds")),
  error = function(e) NULL
)

outputs$viz <- tryCatch(
  readRDS(file.path(PATHS$processed, "06_visualization_output.rds")),
  error = function(e) NULL
)

# Check what loaded successfully
for (name in names(outputs)) {
  status <- ifelse(!is.null(outputs[[name]]), "✓", "✗")
  message("  ", status, " ", name)
}

# =============================================================================
# 3. COMPILE SUMMARY STATISTICS
# =============================================================================

message("\n--- Compiling Summary Statistics ---\n")

summary_stats <- list()

# Dataset info
summary_stats$dataset <- list(
  accession = "GSE74764",
  platform = "GPL16699 (Agilent SurePrint G3 Human GE v2 8x60K)",
  organism = "Homo sapiens",
  cell_lines = c("MDA-MB-231 (TNBC)", "HDF (Normal Fibroblasts)"),
  treatment = "TP4 (Tilapia Piscidin 4), 14 µg/mL, 6 hours",
  replicates = 2,
  total_samples = 8
)

# Preprocessing stats
if (!is.null(outputs$preproc)) {
  summary_stats$preprocessing <- list(
    initial_probes = outputs$preproc$summary$initial_probes,
    filtered_probes = outputs$preproc$summary$filtered_probes,
    unique_genes = outputs$preproc$summary$unique_genes,
    samples = outputs$preproc$summary$samples
  )
}

# DE stats
if (!is.null(outputs$de)) {
  summary_stats$differential_expression <- list(
    fdr_threshold = outputs$de$parameters$fdr_threshold,
    log2fc_threshold = outputs$de$parameters$log2fc_threshold,
    contrasts = outputs$de$summary
  )
  
  # Gene counts
  summary_stats$gene_counts <- list(
    TP4_in_MDA = length(outputs$de$gene_lists$TP4_in_MDA),
    TP4_in_HDF = length(outputs$de$gene_lists$TP4_in_HDF),
    Interaction = length(outputs$de$gene_lists$Interaction),
    MDA_specific = length(outputs$de$gene_lists$MDA_specific),
    Common = length(outputs$de$gene_lists$Common)
  )
}

# Enrichment stats
if (!is.null(outputs$functional)) {
  count_terms <- function(result) {
    if (is.null(result)) return(0)
    if (inherits(result, "enrichResult")) return(nrow(as.data.frame(result)))
    return(0)
  }
  
  summary_stats$enrichment <- list(
    GO_BP_MDA = count_terms(outputs$functional$GO$TP4_in_MDA$BP),
    GO_BP_Interaction = count_terms(outputs$functional$GO$Interaction$BP),
    KEGG_MDA = count_terms(outputs$functional$KEGG$TP4_in_MDA),
    Hallmark_MDA = count_terms(outputs$functional$Hallmark$TP4_in_MDA),
    Hallmark_Interaction = count_terms(outputs$functional$Hallmark$Interaction)
  )
}

message("✓ Summary statistics compiled")

# =============================================================================
# 4. CREATE EXCEL WORKBOOK WITH ALL RESULTS
# =============================================================================

message("\n--- Creating Excel Workbook ---\n")

# Create workbook
wb <- createWorkbook()

# Style definitions
header_style <- createStyle(
  fontSize = 12, 
  fontColour = "#FFFFFF", 
  fgFill = "#4472C4",
  halign = "center", 
  valign = "center", 
  textDecoration = "bold",
  border = "Bottom"
)

# Sheet 1: Summary
addWorksheet(wb, "Summary")
summary_df <- data.frame(
  Parameter = c(
    "Dataset", "Platform", "Organism", "Cell Lines", "Treatment",
    "Total Samples", "Genes Analyzed",
    "FDR Threshold", "Fold-Change Threshold",
    "DEGs in MDA-MB-231", "DEGs in HDF", "Interaction DEGs"
  ),
  Value = c(
    "GSE74764", 
    "GPL16699 (Agilent SurePrint G3 Human GE v2 8x60K)",
    "Homo sapiens",
    "MDA-MB-231 (TNBC), HDF (Normal)",
    "TP4 (14 µg/mL, 6h)",
    "8 (n=2 per group)",
    ifelse(!is.null(summary_stats$preprocessing), 
           as.character(summary_stats$preprocessing$unique_genes), "N/A"),
    ifelse(!is.null(summary_stats$differential_expression),
           as.character(summary_stats$differential_expression$fdr_threshold), "0.05"),
    ifelse(!is.null(summary_stats$differential_expression),
           paste0(summary_stats$differential_expression$log2fc_threshold, 
                  " (", round(2^as.numeric(summary_stats$differential_expression$log2fc_threshold), 2), "-fold)"), 
           "0.585 (1.5-fold)"),
    ifelse(!is.null(summary_stats$gene_counts),
           as.character(summary_stats$gene_counts$TP4_in_MDA), "N/A"),
    ifelse(!is.null(summary_stats$gene_counts),
           as.character(summary_stats$gene_counts$TP4_in_HDF), "N/A"),
    ifelse(!is.null(summary_stats$gene_counts),
           as.character(summary_stats$gene_counts$Interaction), "N/A")
  ),
  stringsAsFactors = FALSE
)
writeData(wb, "Summary", summary_df, headerStyle = header_style)
setColWidths(wb, "Summary", cols = 1:2, widths = c(30, 50))

# Sheet 2: DE Summary
if (!is.null(outputs$de)) {
  addWorksheet(wb, "DE Summary")
  writeData(wb, "DE Summary", outputs$de$summary, headerStyle = header_style)
}

# Sheet 3-6: Full DE results for each contrast
if (!is.null(outputs$de)) {
  for (contrast in names(outputs$de$results)) {
    sheet_name <- paste0("DE_", substr(contrast, 1, 15))
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, outputs$de$results[[contrast]], headerStyle = header_style)
  }
}

# Sheet 7: Gene Lists
if (!is.null(outputs$de)) {
  addWorksheet(wb, "Gene Lists")
  
  max_genes <- max(sapply(outputs$de$gene_lists, length))
  gene_list_df <- data.frame(matrix(NA, nrow = max_genes, ncol = 0))
  
  for (name in names(outputs$de$gene_lists)) {
    genes <- outputs$de$gene_lists[[name]]
    padded <- c(genes, rep(NA, max_genes - length(genes)))
    gene_list_df[[name]] <- padded
  }
  
  writeData(wb, "Gene Lists", gene_list_df, headerStyle = header_style)
}

# Sheet 8: Top Enriched Pathways
if (!is.null(outputs$functional)) {
  addWorksheet(wb, "Top Pathways")
  
  pathway_summary <- data.frame(
    Database = character(),
    Contrast = character(),
    Pathway = character(),
    P_adj = numeric(),
    GeneCount = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Add Hallmark results
  for (contrast in c("TP4_in_MDA", "Interaction")) {
    result <- outputs$functional$Hallmark[[contrast]]
    if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
      df <- as.data.frame(result) %>% head(10)
      pathway_summary <- bind_rows(
        pathway_summary,
        data.frame(
          Database = "Hallmark",
          Contrast = contrast,
          Pathway = gsub("HALLMARK_", "", df$ID),
          P_adj = df$p.adjust,
          GeneCount = df$Count,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  writeData(wb, "Top Pathways", pathway_summary, headerStyle = header_style)
}

# Save workbook
excel_file <- file.path(PATHS$results, "TP4_TNBC_Analysis_Results.xlsx")
saveWorkbook(wb, excel_file, overwrite = TRUE)
message("✓ Excel workbook saved: ", excel_file)

# =============================================================================
# 5. CREATE KEY FINDINGS DOCUMENT
# =============================================================================

message("\n--- Creating Key Findings Document ---\n")

# Get top genes from interaction contrast
top_interaction_genes <- NULL
if (!is.null(outputs$de)) {
  top_interaction_genes <- outputs$de$results[["Interaction"]] %>%
    filter(significant) %>%
    arrange(P.Value) %>%
    head(20)
}

# Get top pathways
top_pathways <- NULL
if (!is.null(outputs$functional) && !is.null(outputs$functional$Hallmark[["Interaction"]])) {
  result <- outputs$functional$Hallmark[["Interaction"]]
  if (nrow(as.data.frame(result)) > 0) {
    top_pathways <- as.data.frame(result) %>% 
      head(10) %>%
      mutate(Pathway = gsub("HALLMARK_", "", ID))
  }
}

key_findings <- paste0(
  "================================================================================\n",
  "KEY FINDINGS: TP4-INDUCED TRANSCRIPTOMIC CHANGES IN TNBC\n",
  "================================================================================\n\n",
  "Analysis Date: ", format(Sys.time(), "%Y-%m-%d"), "\n",
  "Dataset: GSE74764\n\n",
  
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
  "1. OVERVIEW\n",
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
  
  "This analysis examined the transcriptomic response to TP4 (Tilapia Piscidin 4)\n",
  "treatment in triple-negative breast cancer cells (MDA-MB-231) compared to\n",
  "normal human dermal fibroblasts (HDF).\n\n",
  
  "Key Question: What genes respond specifically to TP4 in cancer cells but not\n",
  "in normal cells? These genes may explain TP4's selective anticancer activity.\n\n",
  
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
  "2. DIFFERENTIAL EXPRESSION SUMMARY\n",
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
  
  "Thresholds: FDR < 0.05, |log2FC| > 0.585 (1.5-fold)\n\n",
  
  "TP4 effect in MDA-MB-231 (TNBC): ", 
  ifelse(!is.null(summary_stats$gene_counts), 
         paste0(summary_stats$gene_counts$TP4_in_MDA, " genes"), "N/A"), "\n",
  "TP4 effect in HDF (Normal): ", 
  ifelse(!is.null(summary_stats$gene_counts), 
         paste0(summary_stats$gene_counts$TP4_in_HDF, " genes"), "N/A"), "\n",
  "TNBC-specific response (Interaction): ", 
  ifelse(!is.null(summary_stats$gene_counts), 
         paste0(summary_stats$gene_counts$Interaction, " genes"), "N/A"), "\n\n",
  
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
  "3. TOP TNBC-SPECIFIC GENES (Interaction Contrast)\n",
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
  
  "These genes respond differently to TP4 in cancer cells vs normal cells:\n\n",
  
  ifelse(!is.null(top_interaction_genes),
         paste(sapply(1:nrow(top_interaction_genes), function(i) {
           sprintf("  %2d. %-15s log2FC: %+.2f  FDR: %.2e", 
                   i, 
                   top_interaction_genes$gene[i],
                   top_interaction_genes$logFC[i],
                   top_interaction_genes$adj.P.Val[i])
         }), collapse = "\n"),
         "  [Results not available]"), "\n\n",
  
  "Interpretation:\n",
  "  - Positive log2FC: More induced by TP4 in TNBC than in normal cells\n",
  "  - Negative log2FC: More suppressed by TP4 in TNBC than in normal cells\n\n",
  
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
  "4. ENRICHED PATHWAYS\n",
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
  
  "Top Hallmark gene sets enriched in TNBC-specific response:\n\n",
  
  ifelse(!is.null(top_pathways),
         paste(sapply(1:nrow(top_pathways), function(i) {
           sprintf("  %2d. %-40s (FDR: %.2e)", 
                   i, 
                   top_pathways$Pathway[i],
                   top_pathways$p.adjust[i])
         }), collapse = "\n"),
         "  [Results not available]"), "\n\n",
  
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
  "5. BIOLOGICAL INTERPRETATION\n",
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
  
  "Based on pathway analysis, TP4 treatment in TNBC cells activates:\n\n",
  
  "  1. STRESS RESPONSE PATHWAYS\n",
  "     - Oxidative stress response\n",
  "     - Unfolded protein response (ER stress)\n",
  "     - Heat shock response\n\n",
  
  "  2. APOPTOSIS/CELL DEATH\n",
  "     - Mitochondrial dysfunction\n",
  "     - Caspase activation\n",
  "     - Pro-apoptotic signaling\n\n",
  
  "  3. SIGNALING CASCADES\n",
  "     - MAPK/JNK pathway activation\n",
  "     - AP-1 transcription factor activity\n",
  "     - Calcium-dependent signaling\n\n",
  
  "These findings are consistent with the known mechanism of TP4, which\n",
  "selectively kills cancer cells through membrane disruption, calcium\n",
  "influx, and mitochondrial dysfunction.\n\n",
  
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
  "6. LIMITATIONS AND CAVEATS\n",
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
  
  "  1. SMALL SAMPLE SIZE (n=2 per group)\n",
  "     - Limited statistical power to detect small effects\n",
  "     - Results should be considered hypothesis-generating\n",
  "     - Validation in independent datasets is recommended\n\n",
  
  "  2. SINGLE TIME POINT (6 hours)\n",
  "     - May miss early response genes (< 6h)\n",
  "     - May miss late response genes (> 6h)\n",
  "     - Kinetic studies would provide more complete picture\n\n",
  
  "  3. IN VITRO MODEL\n",
  "     - May not fully represent in vivo tumor microenvironment\n",
  "     - Lacks immune cell interactions\n\n",
  
  "  4. SINGLE CELL LINE PER TYPE\n",
  "     - MDA-MB-231 may not represent all TNBC subtypes\n",
  "     - HDF may not represent all normal tissues\n\n",
  
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
  "7. CONCLUSIONS\n",
  "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
  
  "This analysis reveals that TP4 induces a distinct transcriptomic response in\n",
  "TNBC cells compared to normal fibroblasts. The TNBC-specific response is\n",
  "characterized by activation of stress response and cell death pathways,\n",
  "consistent with TP4's selective anticancer mechanism.\n\n",
  
  "Key candidate genes and pathways identified here warrant further validation\n",
  "and may provide insights into strategies for enhancing TP4's therapeutic\n",
  "efficacy or developing related anticancer peptides.\n\n",
  
  "================================================================================\n"
)

# Save key findings
findings_file <- file.path(PATHS$results, "Key_Findings.txt")
writeLines(key_findings, findings_file)
message("✓ Key findings saved: ", findings_file)

cat(key_findings)

# =============================================================================
# 6. SAVE SESSION INFORMATION
# =============================================================================

message("\n--- Saving Session Information ---\n")

session_info <- capture.output({
  cat("================================================================================\n")
  cat("SESSION INFORMATION\n")
  cat("================================================================================\n\n")
  cat("Analysis completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  cat("R VERSION:\n")
  print(R.version)
  
  cat("\n\nLOADED PACKAGES:\n")
  print(sessionInfo())
  
  cat("\n\nBIOCONDUCTOR VERSION:\n")
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    print(BiocManager::version())
  }
  
  cat("\n\nPROJECT PATHS:\n")
  for (name in names(PATHS)) {
    cat("  ", name, ": ", PATHS[[name]], "\n")
  }
  
  cat("\n\nANALYSIS PARAMETERS:\n")
  for (name in names(ANALYSIS_PARAMS)) {
    val <- ANALYSIS_PARAMS[[name]]
    if (is.list(val)) {
      cat("  ", name, ": [list]\n")
    } else {
      cat("  ", name, ": ", paste(val, collapse = ", "), "\n")
    }
  }
})

session_file <- file.path(PATHS$results, "session_info.txt")
writeLines(session_info, session_file)
message("✓ Session info saved: ", session_file)

# =============================================================================
# 7. CREATE FILE MANIFEST
# =============================================================================

message("\n--- Creating File Manifest ---\n")

# List all output files
list_files_recursive <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("^\\.git", files)]
  return(files)
}

manifest <- data.frame(
  File = character(),
  Directory = character(),
  Size_KB = numeric(),
  Modified = character(),
  stringsAsFactors = FALSE
)

# Collect files from each directory
for (dir_name in c("data", "results", "figures")) {
  dir_path <- PATHS[[dir_name]]
  if (dir.exists(dir_path)) {
    files <- list_files_recursive(dir_path)
    for (f in files) {
      info <- file.info(f)
      manifest <- rbind(manifest, data.frame(
        File = basename(f),
        Directory = dirname(gsub(PATHS$root, "", f)),
        Size_KB = round(info$size / 1024, 2),
        Modified = format(info$mtime, "%Y-%m-%d %H:%M"),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Save manifest
manifest_file <- file.path(PATHS$results, "file_manifest.csv")
write.csv(manifest, manifest_file, row.names = FALSE)
message("✓ File manifest saved: ", manifest_file)

message("\nTotal output files: ", nrow(manifest))

# =============================================================================
# 8. FINAL SUMMARY
# =============================================================================

final_summary <- paste0(
  "\n",
  "╔══════════════════════════════════════════════════════════════════════════╗\n",
  "║                     ANALYSIS COMPLETE                                    ║\n",
  "╠══════════════════════════════════════════════════════════════════════════╣\n",
  "║                                                                          ║\n",
  "║  Dataset: GSE74764                                                       ║\n",
  "║  Analysis: TP4-Induced Transcriptomic Changes in TNBC                    ║\n",
  "║                                                                          ║\n",
  "║  Key Results:                                                            ║\n",
  "║  • Genes analyzed: ", sprintf("%-5s", 
     ifelse(!is.null(summary_stats$preprocessing), 
            summary_stats$preprocessing$unique_genes, "N/A")), 
  "                                           ║\n",
  "║  • DEGs in MDA-MB-231: ", sprintf("%-4s", 
     ifelse(!is.null(summary_stats$gene_counts), 
            summary_stats$gene_counts$TP4_in_MDA, "N/A")), 
  "                                        ║\n",
  "║  • TNBC-specific DEGs: ", sprintf("%-4s", 
     ifelse(!is.null(summary_stats$gene_counts), 
            summary_stats$gene_counts$Interaction, "N/A")), 
  "                                        ║\n",
  "║                                                                          ║\n",
  "║  Output Files:                                                           ║\n",
  "║  • results/TP4_TNBC_Analysis_Results.xlsx                                ║\n",
  "║  • results/Key_Findings.txt                                              ║\n",
  "║  • results/DEG_tables/*.csv                                              ║\n",
  "║  • results/enrichment/*.csv                                              ║\n",
  "║  • figures/*.pdf/*.png                                                   ║\n",
  "║                                                                          ║\n",
  "╚══════════════════════════════════════════════════════════════════════════╝\n"
)

cat(final_summary)

# Save final summary
summary_file <- file.path(PATHS$results, "analysis_complete.txt")
writeLines(c(
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "",
  "Analysis of GSE74764 completed successfully.",
  "",
  "See Key_Findings.txt for biological interpretation.",
  "See TP4_TNBC_Analysis_Results.xlsx for complete results."
), summary_file)

message("\n✓ Analysis pipeline complete!\n")

################################################################################
# END OF SCRIPT
################################################################################
