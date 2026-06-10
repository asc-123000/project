#!/usr/bin/env Rscript
# PDF to PNG Conversion Script for TP4-TNBC Analysis Figures
# This script converts PDF figures to high-quality PNG format
# Requirements: magick or pdftools R packages

cat("\n=== PDF to PNG Conversion for Research Figures ===\n\n")

# Install required packages if needed
required_packages <- c("magick")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing package: %s\n", pkg))
    install.packages(pkg, repos = "http://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# Set base directory
base_dir <- "outputs/Figures"
pdf_files_to_convert <- c(
  "02_Supplementary_Figures/FigS1_DEG_Summary/FigS1_DEG_Summary.pdf",
  "03_QualityControl/BoxPlots/boxplot_by_group.pdf",
  "03_QualityControl/DensityPlots/density_by_cellline.pdf",
  "03_QualityControl/DensityPlots/density_by_treatment.pdf",
  "03_QualityControl/PCA_Plots/PCA_by_group.pdf",
  "03_QualityControl/PCA_Plots/PCA_PC3_PC4.pdf",
  "03_QualityControl/PCA_Plots/scree_plot.pdf",
  "03_QualityControl/SampleCorrelation/sample_dendrogram.pdf",
  "03_QualityControl/VariableGenesHeatmap/heatmap_top_variable_genes.pdf",
  "04_DifferentialExpression/DiagnosticPlots/mean_variance_trend.pdf",
  "04_DifferentialExpression/DiagnosticPlots/pvalue_distributions.pdf",
  "04_DifferentialExpression/MA_Plots/MA_plots.pdf",
  "05_PathwayAnalysis/KEGG_Pathways/KEGG_Interaction.pdf"
)

# Conversion function with quality settings
convert_pdf_to_png <- function(pdf_path, output_path = NULL, density = 300, quality = 95) {
  
  # Check if PDF exists
  if (!file.exists(pdf_path)) {
    warning(sprintf("File not found: %s", pdf_path))
    return(FALSE)
  }
  
  # Set output path if not specified
  if (is.null(output_path)) {
    output_path <- sub("\\.pdf$", ".png", pdf_path)
  }
  
  # Check if PNG already exists
  if (file.exists(output_path)) {
    cat(sprintf("✓ Already exists: %s\n", basename(output_path)))
    return(TRUE)
  }
  
  tryCatch({
    # Read PDF with high quality settings
    cat(sprintf("Converting: %s\n", basename(pdf_path)))
    
    img <- magick::image_read_pdf(pdf_path, density = density)
    
    # Write PNG with compression
    magick::image_write(img, path = output_path, format = "png", quality = quality, compression = "Zip")
    
    cat(sprintf("✓ Success: %s\n", basename(output_path)))
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("✗ Error converting %s: %s\n", basename(pdf_path), e$message))
    return(FALSE)
  })
}

# Execute conversions
cat("Starting PDF to PNG conversions...\n")
cat("Resolution: 300 DPI | Quality: 95%\n\n")

success_count <- 0
total_count <- length(pdf_files_to_convert)

for (pdf_file in pdf_files_to_convert) {
  full_path <- file.path(base_dir, pdf_file)
  if (convert_pdf_to_png(full_path, density = 300, quality = 95)) {
    success_count <- success_count + 1
  }
}

# Summary report
cat("\n=== Conversion Summary ===\n")
cat(sprintf("Successfully converted: %d/%d files\n", success_count, total_count))
cat(sprintf("Completion rate: %.1f%%\n\n", (success_count/total_count)*100)

# List PNG files created
cat("PNG files available in outputs/Figures:\n")
png_files <- list.files(base_dir, pattern = "\\.png$", recursive = TRUE)
cat(sprintf("Total PNG files: %d\n\n", length(png_files)))

cat("✓ Conversion process complete!\n\n")
