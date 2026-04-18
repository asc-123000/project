################################################################################
# run_analysis.R - Master Script to Run Complete Analysis Pipeline
################################################################################
#
# Project: TP4-Induced Transcriptomic Changes in Triple-Negative Breast Cancer
# 
# INSTRUCTIONS FOR WINDOWS / RSTUDIO:
# ─────────────────────────────────────────────────
#   1. Open this file in RStudio
#   2. Click "Source" button (top-right of script editor), OR
#   3. Run: source("run_analysis.R") from the console
#
# PREREQUISITES:
#   - R version 4.0+ recommended
#   - Internet connection for downloading packages and GEO data
#   - On Windows: Rtools may be required for some Bioconductor packages
#     Download from: https://cran.r-project.org/bin/windows/Rtools/
#
# Author: Computational Biology Pipeline
# Date: 2026-01-19
#
################################################################################

# =============================================================================
# 1. SET UP WORKING DIRECTORY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("  TP4 TNBC Transcriptomics Analysis Pipeline - Windows/RStudio Compatible\n")
cat("================================================================================\n")
cat("\n")

# Determine project root and set working directory
find_project_root <- function() {
  # Method 1: Use 'here' package if available
  if (requireNamespace("here", quietly = TRUE)) {
    root <- here::here()
    if (file.exists(file.path(root, "scripts", "00_setup.R"))) {
      return(normalizePath(root, winslash = "/"))
    }
  }
  
  # Method 2: Check if running in RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    script_path <- tryCatch({
      rstudioapi::getActiveDocumentContext()$path
    }, error = function(e) "")
    
    if (nzchar(script_path)) {
      dir <- dirname(script_path)
      if (file.exists(file.path(dir, "scripts", "00_setup.R"))) {
        return(normalizePath(dir, winslash = "/"))
      }
    }
  }
  
  # Method 3: Current directory
  cwd <- getwd()
  if (file.exists(file.path(cwd, "scripts", "00_setup.R"))) {
    return(normalizePath(cwd, winslash = "/"))
  }
  
  stop("Cannot find project root. Please open this script from the project folder.")
}

project_root <- find_project_root()
setwd(project_root)
cat("Working directory set to:", project_root, "\n\n")

# =============================================================================
# 2. CHECK WINDOWS REQUIREMENTS
# =============================================================================

check_windows_requirements <- function() {
  if (.Platform$OS.type == "windows") {
    cat("Detected Windows operating system.\n")
    
    # Check for Rtools
    rtools_installed <- FALSE
    
    # Check via pkgbuild if available
    if (requireNamespace("pkgbuild", quietly = TRUE)) {
      rtools_installed <- pkgbuild::has_rtools()
    } else {
      # Manual check for common Rtools paths
      rtools_paths <- c(
        "C:/rtools44",
        "C:/rtools43",
        "C:/rtools42",
        "C:/rtools40",
        "C:/Rtools"
      )
      rtools_installed <- any(sapply(rtools_paths, dir.exists))
    }
    
    if (!rtools_installed) {
      cat("\n")
      cat("WARNING: Rtools may not be installed!\n")
      cat("Some Bioconductor packages require Rtools for compilation.\n")
      cat("Download from: https://cran.r-project.org/bin/windows/Rtools/\n")
      cat("\n")
      
      answer <- readline("Continue anyway? (y/n): ")
      if (tolower(answer) != "y") {
        stop("Aborted. Please install Rtools and try again.")
      }
    } else {
      cat("Rtools detected.\n")
    }
  }
}

check_windows_requirements()

# =============================================================================
# 3. INSTALL/LOAD REQUIRED PACKAGES
# =============================================================================

cat("\n--- Setting up environment ---\n\n")

# Source the setup script
source("scripts/00_setup.R")

# Check if packages need installation
all_packages <- c(cran_packages, bioc_packages)
missing <- all_packages[!(all_packages %in% installed.packages()[, "Package"])]

if (length(missing) > 0) {
  cat("\nMissing", length(missing), "packages. Installing...\n")
  cat("This may take 10-30 minutes on first run.\n\n")
  
  install_all_packages()
}

# Load all packages
cat("\nLoading packages...\n")
load_all_packages(verbose = FALSE)
cat("All packages loaded.\n")

# =============================================================================
# 4. RUN ANALYSIS PIPELINE
# =============================================================================

run_pipeline <- function(scripts = NULL) {
  # Default: run all scripts in order
  if (is.null(scripts)) {
    scripts <- c(
      "01_data_download.R",
      "02_preprocessing.R",
      "03_quality_control.R",
      "04_differential_expression.R",
      "05_functional_analysis.R",
      "06_visualization.R",
      "07_generate_report.R"
    )
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("  RUNNING ANALYSIS PIPELINE\n")
  cat("================================================================================\n")
  cat("\n")
  
  for (script in scripts) {
    script_path <- file.path("scripts", script)
    
    if (!file.exists(script_path)) {
      cat("ERROR: Script not found:", script_path, "\n")
      next
    }
    
    cat("\n")
    cat("--------------------------------------------------------------------------------\n")
    cat("Running:", script, "\n")
    cat("--------------------------------------------------------------------------------\n")
    cat("\n")
    
    # Run script with error handling
    result <- tryCatch({
      source(script_path, local = new.env())
      "SUCCESS"
    }, error = function(e) {
      cat("\nERROR in", script, ":\n")
      cat(conditionMessage(e), "\n")
      return("FAILED")
    }, warning = function(w) {
      cat("WARNING:", conditionMessage(w), "\n")
      return("SUCCESS with warnings")
    })
    
    cat("\n", script, ":", result, "\n")
    
    if (result == "FAILED") {
      answer <- readline("Continue with next script? (y/n): ")
      if (tolower(answer) != "y") {
        cat("Pipeline stopped.\n")
        return(invisible(FALSE))
      }
    }
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("  PIPELINE COMPLETE\n")
  cat("================================================================================\n")
  cat("\n")
  cat("Results saved in:\n")
  cat("  - Data:    ", file.path(project_root, "data", "processed"), "\n")
  cat("  - Figures: ", file.path(project_root, "figures"), "\n")
  cat("  - Tables:  ", file.path(project_root, "results"), "\n")
  cat("\n")
  
  return(invisible(TRUE))
}

# =============================================================================
# 5. INTERACTIVE MENU
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("  READY TO RUN\n")
cat("================================================================================\n")
cat("\n")
cat("Options:\n")
cat("  1. run_pipeline()           - Run complete analysis (all scripts)\n")
cat("  2. run_pipeline(c('01_data_download.R'))  - Run specific script(s)\n")
cat("  3. Source individual scripts from the 'scripts/' folder\n")
cat("\n")
cat("To start the full analysis, run:\n")
cat("  run_pipeline()\n")
cat("\n")

################################################################################
# END OF SCRIPT
################################################################################
