################################################################################
# 00_setup.R - Environment Setup and Package Installation
################################################################################
#
# Project: TP4-Induced Transcriptomic Changes in Triple-Negative Breast Cancer
# Dataset: GSE74764 (Agilent SurePrint G3 Human GE v2 8x60K)
# 
# Purpose:
#   This script installs and loads all required packages for the complete
#   analysis pipeline. It also sets up global options and defines project paths.
#
# Scientific Rationale:
#   Reproducible research requires explicit documentation of the computational
#   environment. This script ensures all dependencies are available and 
#   establishes consistent settings across all analysis scripts.
#
# Author: Computational Biology Pipeline
# Date: 2026-01-19
#
################################################################################

# =============================================================================
# 1. SET GLOBAL OPTIONS
# =============================================================================

# Set CRAN mirror for reproducibility
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Prevent scientific notation for better readability
options(scipen = 999)

# Set seed for reproducibility (important for any stochastic processes)
set.seed(42)

# =============================================================================
# 2. DEFINE PROJECT PATHS
# =============================================================================

# Get the script directory and set project root
# This allows the pipeline to run from any working directory
get_project_root <- function() {
  # Try to get script path if running via source()
  script_path <- tryCatch({
    dirname(sys.frame(1)$ofile)
  }, error = function(e) {
    # If that fails, use current working directory
    getwd()
  })
  
  # Navigate up from scripts/ to project root
  if (basename(script_path) == "scripts") {
    return(dirname(script_path))
  } else {
    # Assume we're in project root
    return(script_path)
  }
}

# Define project structure
setup_project_paths <- function(project_root = NULL) {
  if (is.null(project_root)) {
    project_root <- get_project_root()
  }
  
  paths <- list(
    root = project_root,
    data = file.path(project_root, "data"),
    raw = file.path(project_root, "data", "raw"),
    processed = file.path(project_root, "data", "processed"),
    scripts = file.path(project_root, "scripts"),
    results = file.path(project_root, "results"),
    deg_tables = file.path(project_root, "results", "DEG_tables"),
    enrichment = file.path(project_root, "results", "enrichment"),
    figures = file.path(project_root, "figures"),
    figures_qc = file.path(project_root, "figures", "QC"),
    figures_de = file.path(project_root, "figures", "DE"),
    figures_pathways = file.path(project_root, "figures", "pathways"),
    docs = file.path(project_root, "docs")
  )
  
  # Create directories if they don't exist
  lapply(paths, function(p) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE)
  })
  
  return(paths)
}

# =============================================================================
# 3. PACKAGE INSTALLATION FUNCTIONS
# =============================================================================

#' Install CRAN packages if not already installed
#' @param packages Character vector of package names
install_cran_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    message("Installing CRAN packages: ", paste(new_packages, collapse = ", "))
    install.packages(new_packages, dependencies = TRUE)
  } else {
    message("All CRAN packages already installed.")
  }
}

#' Install Bioconductor packages if not already installed
#' @param packages Character vector of package names
install_bioc_packages <- function(packages) {
  # Ensure BiocManager is available
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    message("Installing Bioconductor packages: ", paste(new_packages, collapse = ", "))
    BiocManager::install(new_packages, ask = FALSE, update = FALSE)
  } else {
    message("All Bioconductor packages already installed.")
  }
}

# =============================================================================
# 4. REQUIRED PACKAGES
# =============================================================================

# CRAN packages
cran_packages <- c(
  # Data manipulation
  "tidyverse",      # Core tidyverse (dplyr, ggplot2, tidyr, readr, etc.)
  "data.table",     # Fast data manipulation
  
  # Visualization
  "pheatmap",       # Heatmaps
  "RColorBrewer",   # Color palettes
  "viridis",        # Color-blind friendly palettes
  "ggrepel",        # Non-overlapping text labels
  "gridExtra",      # Arrange multiple plots
  "scales",         # Scale functions for ggplot2
  "cowplot",        # Publication-ready plots
  
  # Statistical
  "broom",          # Tidy statistical results
  "ggpubr",         # Publication-ready ggplots
  
  # Utilities
  "here",           # Project-relative paths
  "openxlsx",       # Excel file handling
  "knitr",          # Report generation
  "rmarkdown",      # Markdown reports
  "DT"              # Interactive tables
)

# Bioconductor packages
bioc_packages <- c(
  # Data access
  "GEOquery",       # Download GEO datasets
  
  # Microarray analysis
  "limma",          # Linear models for microarray
  "affy",           # Affymetrix analysis (general microarray utilities)
  
  # Annotation
  "AnnotationDbi",  # Annotation database interface
  "org.Hs.eg.db",   # Human gene annotations
  "hgu133plus2.db", # Agilent annotation (fallback)
  
  # Functional analysis
  "clusterProfiler", # GO/KEGG enrichment
  "enrichplot",      # Visualization for enrichment
  "DOSE",            # Disease ontology
  "ReactomePA",      # Reactome pathway analysis
  "msigdbr",         # MSigDB gene sets
  "fgsea",           # Fast GSEA
  "GO.db",           # GO database
  
  # Quality control
  "arrayQualityMetrics", # Array QC (optional, may have dependencies)
  
  # Utilities
  "Biobase"          # Core Bioconductor classes
)

# =============================================================================
# 5. INSTALL ALL PACKAGES
# =============================================================================

install_all_packages <- function() {
  message("\n", paste(rep("=", 60), collapse = ""))
  message("INSTALLING REQUIRED PACKAGES")
  message(paste(rep("=", 60), collapse = ""), "\n")
  
  # Install CRAN packages
  message("\n--- CRAN Packages ---\n")
  install_cran_packages(cran_packages)
  
  # Install Bioconductor packages
  message("\n--- Bioconductor Packages ---\n")
  install_bioc_packages(bioc_packages)
  
  message("\n", paste(rep("=", 60), collapse = ""))
  message("PACKAGE INSTALLATION COMPLETE")
  message(paste(rep("=", 60), collapse = ""), "\n")
}

# =============================================================================
# 6. LOAD ALL PACKAGES
# =============================================================================

load_all_packages <- function(verbose = TRUE) {
  all_packages <- c(cran_packages, bioc_packages)
  
  if (verbose) {
    message("\n--- Loading packages ---\n")
  }
  
  # Suppress startup messages for cleaner output
  loaded <- sapply(all_packages, function(pkg) {
    tryCatch({
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      if (verbose) message("  ✓ ", pkg)
      TRUE
    }, error = function(e) {
      if (verbose) message("  ✗ ", pkg, " - ", e$message)
      FALSE
    })
  })
  
  if (verbose) {
    message("\nSuccessfully loaded: ", sum(loaded), "/", length(all_packages), " packages")
  }
  
  invisible(loaded)
}

# =============================================================================
# 7. SESSION INFO FUNCTION
# =============================================================================

#' Save session information for reproducibility
#' @param output_dir Directory to save session info
save_session_info <- function(output_dir) {
  session_file <- file.path(output_dir, "session_info.txt")
  
  sink(session_file)
  cat("================================================================================\n")
  cat("SESSION INFORMATION\n")
  cat("================================================================================\n\n")
  cat("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("R Version:\n")
  print(R.version)
  cat("\n\nLoaded Packages:\n")
  print(sessionInfo())
  cat("\n\nBioconductor Version:\n")
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    print(BiocManager::version())
  }
  sink()
  
  message("Session info saved to: ", session_file)
}

# =============================================================================
# 8. ANALYSIS PARAMETERS
# =============================================================================

# Define global analysis parameters as a list for consistency across scripts
ANALYSIS_PARAMS <- list(
  # Dataset identifiers
  geo_accession = "GSE74764",
  platform = "GPL16699",
  
  # Statistical thresholds
  # FDR < 0.05: Standard threshold, appropriate for discovery given small n
  fdr_threshold = 0.05,
  
  # |log2FC| > 0.585 corresponds to 1.5-fold change
  # Rationale: With n=2 per group, we have limited power to detect small effects
  # 1.5-fold captures biologically meaningful changes while being less stringent
  # than the typical 2-fold cutoff (log2FC > 1)
  log2fc_threshold = 0.585,

  # Sample metadata
  cell_lines = c("MDA-MB-231", "HDF"),
  treatments = c("Mock", "TP4"),
  
  # Biological focus areas (for pathway analysis prioritization)
  focus_pathways = c(
    "apoptosis",
    "cell death",
    "stress response",
    "MAPK signaling",
    "calcium signaling",
    "mitochondria",
    "oxidative stress",
    "ER stress",
    "autophagy"
  ),
  
  # Visualization settings
  plot_width = 8,
  plot_height = 6,
  plot_dpi = 300,
  
  # Color schemes
  colors = list(
    cell_line = c("MDA-MB-231" = "#E41A1C", "HDF" = "#377EB8"),
    treatment = c("Mock" = "#999999", "TP4" = "#4DAF4A"),
    direction = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "#999999")
  )
)

# =============================================================================
# 9. CUSTOM THEME FOR PUBLICATION-QUALITY PLOTS
# =============================================================================

#' Publication-ready ggplot2 theme
#' Based on theme_classic with modifications for journal figures
theme_publication <- function(base_size = 12, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # Title and labels
      plot.title = element_text(size = base_size * 1.2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size * 0.9, color = "black"),
      
      # Legend
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size * 0.9),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      
      # Panel
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Axis lines
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      
      # Margins
      plot.margin = margin(10, 10, 10, 10)
    )
}

# =============================================================================
# 10. MAIN EXECUTION
# =============================================================================

# This section runs when the script is sourced
if (interactive() || !exists(".setup_complete")) {
  message("\n")
  message("╔══════════════════════════════════════════════════════════════════════════╗")
  message("║  TP4 TNBC Transcriptomics Analysis Pipeline - Environment Setup          ║")
  message("╚══════════════════════════════════════════════════════════════════════════╝")
  message("\n")
  
  # Set up project paths
  PATHS <- setup_project_paths()
  message("Project root: ", PATHS$root)
  
  # Check if packages need to be installed
  message("\nChecking package availability...")
  
  all_packages <- c(cran_packages, bioc_packages)
  missing <- all_packages[!(all_packages %in% installed.packages()[, "Package"])]
  
  if (length(missing) > 0) {
    message("\nMissing packages detected: ", length(missing))
    message("Run install_all_packages() to install them.\n")
    message("Missing: ", paste(head(missing, 10), collapse = ", "), 
            if(length(missing) > 10) paste0(" ... and ", length(missing) - 10, " more"))
  } else {
    message("All required packages are installed.")
    message("Loading packages...")
    load_all_packages(verbose = FALSE)
    message("All packages loaded successfully.")
  }
  
  # Mark setup as complete
  .setup_complete <- TRUE
  
  message("\n--- Setup complete ---\n")
  message("Available objects:")
  message("  PATHS           - Project directory paths")
  message("  ANALYSIS_PARAMS - Analysis parameters and thresholds")
  message("  theme_publication() - ggplot2 theme for figures")
  message("\nAvailable functions:")
  message("  install_all_packages() - Install all required packages")
  message("  load_all_packages()    - Load all packages")
  message("  save_session_info()    - Save session info for reproducibility")
  message("\n")
}

################################################################################
# END OF SCRIPT
################################################################################
