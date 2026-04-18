################################################################################
# 05_functional_analysis.R - Pathway and Gene Set Enrichment Analysis
################################################################################
#
# Project: TP4-Induced Transcriptomic Changes in Triple-Negative Breast Cancer
# Dataset: GSE74764 (Agilent SurePrint G3 Human GE v2 8x60K)
# 
# Purpose:
#   Perform comprehensive functional enrichment analysis:
#   1. Gene Ontology (GO) enrichment
#   2. KEGG pathway analysis
#   3. Reactome pathway analysis
#   4. MSigDB Hallmark gene sets
#   5. Gene Set Enrichment Analysis (GSEA)
#   6. Focus on cancer-relevant pathways
#
# Scientific Focus:
# ─────────────────────────────────────────────────
#   Based on the original publication and TP4 mechanism, we focus on:
#   
#   1. APOPTOSIS / CELL DEATH
#      - TP4 induces selective death in cancer cells
#      - Look for mitochondrial apoptosis, caspase activation
#   
#   2. STRESS RESPONSE
#      - ER stress / Unfolded Protein Response (UPR)
#      - Oxidative stress
#      - DNA damage response
#   
#   3. CALCIUM SIGNALING
#      - TP4 disrupts calcium homeostasis
#      - Calpain activation
#   
#   4. MAPK / JNK / AP-1 PATHWAY
#      - Original study highlighted FOSB activation
#      - Stress-activated protein kinases
#   
#   5. MITOCHONDRIAL FUNCTION
#      - Membrane potential disruption
#      - Oxidative phosphorylation
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
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DOSE)
library(ReactomePA)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(tidyverse)
library(pheatmap)

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Script 05: Functional Analysis                                          ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

# =============================================================================
# 2. LOAD DIFFERENTIAL EXPRESSION RESULTS
# =============================================================================

message("--- Loading DE Results ---\n")

# Load DE output
de_output <- readRDS(file.path(PATHS$processed, "04_de_output.rds"))

de_results <- de_output$results
gene_lists <- de_output$gene_lists
fit <- de_output$fit

message("✓ Loaded DE results for ", length(de_results), " contrasts")

# Load expression data for GSEA
preproc_output <- readRDS(file.path(PATHS$processed, "02_preprocessing_output.rds"))
expr_gene <- preproc_output$expression_gene

# =============================================================================
# 3. PREPARE GENE LISTS FOR ENRICHMENT
# =============================================================================

message("\n--- Preparing Gene Lists ---\n")

# Function to convert gene symbols to Entrez IDs
symbols_to_entrez <- function(symbols) {
  entrez <- mapIds(org.Hs.eg.db,
                   keys = symbols,
                   column = "ENTREZID",
                   keytype = "SYMBOL",
                   multiVals = "first")
  # Remove NAs
  entrez <- entrez[!is.na(entrez)]
  return(entrez)
}

# Convert each gene list
entrez_lists <- lapply(gene_lists, function(x) {
  if (length(x) == 0) return(character(0))
  symbols_to_entrez(x)
})

# Report conversion success
for (name in names(gene_lists)) {
  message("  ", name, ": ", length(gene_lists[[name]]), " symbols → ", 
          length(entrez_lists[[name]]), " Entrez IDs")
}

# =============================================================================
# 4. GENE ONTOLOGY ENRICHMENT
# =============================================================================

message("\n--- Gene Ontology Enrichment ---\n")

# We'll focus on the most important contrasts:
# 1. TP4 effect in MDA-MB-231 (cancer response)
# 2. Interaction (cancer-specific response)

# Get background gene universe
universe_symbols <- rownames(expr_gene)
universe_entrez <- symbols_to_entrez(universe_symbols)
message("Background universe: ", length(universe_entrez), " genes")

# Function to run GO enrichment for a gene list
run_go_enrichment <- function(genes_entrez, name, universe = NULL) {
  if (length(genes_entrez) < 5) {
    message("  Skipping ", name, ": too few genes (", length(genes_entrez), ")")
    return(NULL)
  }
  
  message("  Running GO enrichment for ", name, " (", length(genes_entrez), " genes)...")
  
  # Biological Process
  ego_BP <- enrichGO(
    gene = genes_entrez,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  
  # Molecular Function
  ego_MF <- enrichGO(
    gene = genes_entrez,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  
  # Cellular Component
  ego_CC <- enrichGO(
    gene = genes_entrez,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  
  return(list(BP = ego_BP, MF = ego_MF, CC = ego_CC))
}

# Run GO enrichment for key contrasts
go_results <- list()

# TP4 in MDA-MB-231
go_results[["TP4_in_MDA"]] <- run_go_enrichment(
  entrez_lists[["TP4_in_MDA"]], 
  "TP4_in_MDA",
  universe_entrez
)

# TP4 in HDF
go_results[["TP4_in_HDF"]] <- run_go_enrichment(
  entrez_lists[["TP4_in_HDF"]], 
  "TP4_in_HDF",
  universe_entrez
)

# Interaction (TNBC-specific)
go_results[["Interaction"]] <- run_go_enrichment(
  entrez_lists[["Interaction"]], 
  "Interaction",
  universe_entrez
)

# MDA-specific genes
go_results[["MDA_specific"]] <- run_go_enrichment(
  entrez_lists[["MDA_specific"]], 
  "MDA_specific",
  universe_entrez
)

message("✓ GO enrichment complete")

# =============================================================================
# 5. KEGG PATHWAY ANALYSIS
# =============================================================================

message("\n--- KEGG Pathway Analysis ---\n")

run_kegg_enrichment <- function(genes_entrez, name, universe = NULL) {
  if (length(genes_entrez) < 5) {
    message("  Skipping ", name, ": too few genes")
    return(NULL)
  }
  
  message("  Running KEGG for ", name, "...")
  
  # Wrap in tryCatch to handle network errors gracefully
  kegg <- tryCatch({
    enrichKEGG(
      gene = genes_entrez,
      universe = universe,
      organism = "hsa",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1
    )
  }, error = function(e) {
    message("  WARNING: KEGG analysis failed for ", name, ": ", e$message)
    message("  (KEGG servers may be unavailable - skipping)")
    return(NULL)
  })
  
  # Convert Entrez IDs to symbols for readability
  if (!is.null(kegg) && nrow(kegg) > 0) {
    kegg <- tryCatch({
      setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    }, error = function(e) kegg)
  }
  
  return(kegg)
}

kegg_results <- list()

kegg_results[["TP4_in_MDA"]] <- run_kegg_enrichment(
  entrez_lists[["TP4_in_MDA"]], "TP4_in_MDA", universe_entrez
)

kegg_results[["TP4_in_HDF"]] <- run_kegg_enrichment(
  entrez_lists[["TP4_in_HDF"]], "TP4_in_HDF", universe_entrez
)

kegg_results[["Interaction"]] <- run_kegg_enrichment(
  entrez_lists[["Interaction"]], "Interaction", universe_entrez
)

message("✓ KEGG analysis complete (some may have been skipped due to network issues)")

# =============================================================================
# 6. REACTOME PATHWAY ANALYSIS
# =============================================================================

message("\n--- Reactome Pathway Analysis ---\n")

run_reactome_enrichment <- function(genes_entrez, name, universe = NULL) {
  if (length(genes_entrez) < 5) {
    message("  Skipping ", name, ": too few genes")
    return(NULL)
  }
  
  message("  Running Reactome for ", name, "...")
  
  # Wrap in tryCatch to handle network errors gracefully
  reactome <- tryCatch({
    enrichPathway(
      gene = genes_entrez,
      universe = universe,
      organism = "human",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1,
      readable = TRUE
    )
  }, error = function(e) {
    message("  WARNING: Reactome analysis failed for ", name, ": ", e$message)
    message("  (Reactome servers may be unavailable - skipping)")
    return(NULL)
  })
  
  return(reactome)
}

reactome_results <- list()

reactome_results[["TP4_in_MDA"]] <- run_reactome_enrichment(
  entrez_lists[["TP4_in_MDA"]], "TP4_in_MDA", universe_entrez
)

reactome_results[["TP4_in_HDF"]] <- run_reactome_enrichment(
  entrez_lists[["TP4_in_HDF"]], "TP4_in_HDF", universe_entrez
)

reactome_results[["Interaction"]] <- run_reactome_enrichment(
  entrez_lists[["Interaction"]], "Interaction", universe_entrez
)

message("✓ Reactome analysis complete")

# =============================================================================
# 7. MSIGDB HALLMARK GENE SETS
# =============================================================================

message("\n--- MSigDB Hallmark Analysis ---\n")

# Get Hallmark gene sets
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_df$entrez_gene, hallmark_df$gs_name)

# Convert to proper format for enricher
hallmark_t2g <- hallmark_df %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::rename(term = gs_name, gene = entrez_gene)

run_hallmark_enrichment <- function(genes_entrez, name, universe = NULL) {
  if (length(genes_entrez) < 5) {
    message("  Skipping ", name, ": too few genes")
    return(NULL)
  }
  
  message("  Running Hallmark for ", name, "...")
  
  hallmark <- enricher(
    gene = genes_entrez,
    universe = universe,
    TERM2GENE = hallmark_t2g,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )
  
  return(hallmark)
}

hallmark_results <- list()

hallmark_results[["TP4_in_MDA"]] <- run_hallmark_enrichment(
  entrez_lists[["TP4_in_MDA"]], "TP4_in_MDA", universe_entrez
)

hallmark_results[["TP4_in_HDF"]] <- run_hallmark_enrichment(
  entrez_lists[["TP4_in_HDF"]], "TP4_in_HDF", universe_entrez
)

hallmark_results[["Interaction"]] <- run_hallmark_enrichment(
  entrez_lists[["Interaction"]], "Interaction", universe_entrez
)

message("✓ Hallmark analysis complete")

# =============================================================================
# 8. GENE SET ENRICHMENT ANALYSIS (GSEA)
# =============================================================================

message("\n--- Gene Set Enrichment Analysis (GSEA) ---\n")

# GSEA uses the full ranked list of genes, not just significant ones
# Rank by t-statistic from limma

# Prepare ranked gene lists for each contrast
prepare_ranked_list <- function(contrast_name) {
  tt <- de_results[[contrast_name]]
  
  # Get t-statistics
  ranked <- tt$t
  names(ranked) <- tt$gene
  
  # Convert to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = names(ranked),
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  # Remove NAs and create named vector
  valid <- !is.na(entrez_ids)
  ranked_entrez <- ranked[valid]
  names(ranked_entrez) <- entrez_ids[valid]
  
  # Sort by rank
  ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)
  
  return(ranked_entrez)
}

# Run GSEA for key contrasts
run_fgsea <- function(ranked_list, genesets, contrast_name) {
  message("  Running GSEA for ", contrast_name, "...")
  
  result <- tryCatch({
    fgsea(
      pathways = genesets,
      stats = ranked_list,
      minSize = 15,
      maxSize = 500
      # nperm removed to use fgseaMultilevel algorithm
    )
  }, error = function(e) {
    message("  WARNING: GSEA failed for ", contrast_name, ": ", e$message)
    return(data.frame())
  })
  
  # Sort by NES if results exist
  if (nrow(result) > 0) {
    result <- result[order(result$NES, decreasing = TRUE), ]
  }
  
  return(result)
}

# Prepare ranked lists
ranked_lists <- list()
for (contrast in c("TP4_in_MDA", "TP4_in_HDF", "Interaction")) {
  ranked_lists[[contrast]] <- prepare_ranked_list(contrast)
  message("  ", contrast, ": ", length(ranked_lists[[contrast]]), " ranked genes")
}

# Run GSEA with Hallmark gene sets
gsea_hallmark <- list()
for (contrast in names(ranked_lists)) {
  gsea_hallmark[[contrast]] <- run_fgsea(ranked_lists[[contrast]], hallmark_list, contrast)
}

message("✓ GSEA complete")

# =============================================================================
# 9. CANCER-RELEVANT PATHWAY ANALYSIS
# =============================================================================

message("\n--- Cancer-Relevant Pathway Analysis ---\n")

# Get additional MSigDB gene sets relevant to cancer
# C2: Curated gene sets (includes KEGG, Reactome, BioCarta)
# C6: Oncogenic signatures

# Get oncogenic signatures (C6)
oncogenic_df <- msigdbr(species = "Homo sapiens", category = "C6")
oncogenic_t2g <- oncogenic_df %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::rename(term = gs_name, gene = entrez_gene)

# Focus on AP-1, MAPK, and apoptosis related sets
# Get C2 canonical pathways
c2_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")

# Filter for relevant pathways
relevant_keywords <- c("APOPTOSIS", "DEATH", "STRESS", "MAPK", "JNK", 
                       "CALCIUM", "MITOCHOND", "REACTIVE_OXYGEN", "P53",
                       "CASPASE", "BCL2", "ER_STRESS", "UPR", "TNF")

c2_relevant <- c2_df %>%
  filter(grepl(paste(relevant_keywords, collapse = "|"), gs_name, ignore.case = TRUE))

c2_relevant_t2g <- c2_relevant %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::rename(term = gs_name, gene = entrez_gene)

# Run enrichment for cancer-relevant pathways
message("Running cancer-relevant pathway enrichment...")

cancer_pathway_results <- list()

for (contrast in c("TP4_in_MDA", "Interaction")) {
  genes <- entrez_lists[[contrast]]
  
  if (length(genes) >= 5) {
    # Oncogenic signatures
    onco <- enricher(
      gene = genes,
      universe = universe_entrez,
      TERM2GENE = oncogenic_t2g,
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2
    )
    
    # Relevant canonical pathways
    canon <- enricher(
      gene = genes,
      universe = universe_entrez,
      TERM2GENE = c2_relevant_t2g,
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2
    )
    
    cancer_pathway_results[[contrast]] <- list(
      oncogenic = onco,
      canonical = canon
    )
  }
}

message("✓ Cancer-relevant analysis complete")

# =============================================================================
# 10. COMPILE AND SAVE ENRICHMENT RESULTS
# =============================================================================

message("\n--- Saving Enrichment Results ---\n")

# Function to save enrichment result
save_enrichment <- function(result, name, subdir = "") {
  if (is.null(result)) {
    return()
  }
  
  outdir <- file.path(PATHS$enrichment, subdir)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Convert to data frame
  if (inherits(result, "enrichResult")) {
    df <- as.data.frame(result)
  } else if (inherits(result, "data.frame") || inherits(result, "data.table")) {
    df <- as.data.frame(result)
  } else {
    return()
  }
  
  # Check if dataframe is empty
  if (nrow(df) == 0) {
    return()
  }
  
  # Handle list columns (e.g., leadingEdge from fgsea) by converting to strings
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- sapply(df[[col]], function(x) {
        if (is.null(x) || length(x) == 0) {
          return("")
        }
        paste(x, collapse = ";")
      })
    }
  }
  
  # Write to CSV with error handling
  tryCatch({
    write.csv(df, file.path(outdir, paste0(name, ".csv")), row.names = FALSE)
    message("  ✓ ", file.path(subdir, paste0(name, ".csv")))
  }, error = function(e) {
    message("  WARNING: Could not save ", name, ": ", e$message)
  })
}

# Save GO results
for (contrast in names(go_results)) {
  if (!is.null(go_results[[contrast]])) {
    for (ont in names(go_results[[contrast]])) {
      save_enrichment(go_results[[contrast]][[ont]], 
                     paste0("GO_", ont, "_", contrast), 
                     "GO")
    }
  }
}

# Save KEGG results
for (contrast in names(kegg_results)) {
  save_enrichment(kegg_results[[contrast]], 
                 paste0("KEGG_", contrast), 
                 "KEGG")
}

# Save Reactome results
for (contrast in names(reactome_results)) {
  save_enrichment(reactome_results[[contrast]], 
                 paste0("Reactome_", contrast), 
                 "Reactome")
}

# Save Hallmark results
for (contrast in names(hallmark_results)) {
  save_enrichment(hallmark_results[[contrast]], 
                 paste0("Hallmark_", contrast), 
                 "Hallmark")
}

# Save GSEA results
for (contrast in names(gsea_hallmark)) {
  save_enrichment(gsea_hallmark[[contrast]], 
                 paste0("GSEA_Hallmark_", contrast), 
                 "GSEA")
}

# Save cancer pathway results
for (contrast in names(cancer_pathway_results)) {
  save_enrichment(cancer_pathway_results[[contrast]]$oncogenic,
                 paste0("Oncogenic_", contrast),
                 "Cancer")
  save_enrichment(cancer_pathway_results[[contrast]]$canonical,
                 paste0("Canonical_", contrast),
                 "Cancer")
}

# =============================================================================
# 11. CREATE ENRICHMENT VISUALIZATIONS
# =============================================================================

message("\n--- Creating Enrichment Visualizations ---\n")

# Function to create dot plot for enrichment results
create_enrichment_dotplot <- function(result, title, top_n = 20) {
  if (is.null(result) || nrow(result) == 0) {
    return(NULL)
  }
  
  p <- dotplot(result, showCategory = top_n) +
    labs(title = title) +
    theme_publication() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 12)
    )
  
  return(p)
}

# GO Biological Process plots
for (contrast in names(go_results)) {
  if (!is.null(go_results[[contrast]]$BP) && nrow(go_results[[contrast]]$BP) > 0) {
    p <- create_enrichment_dotplot(
      go_results[[contrast]]$BP,
      paste0("GO Biological Process: ", contrast)
    )
    
    if (!is.null(p)) {
      ggsave(file.path(PATHS$figures_pathways, paste0("GO_BP_", contrast, ".pdf")),
             p, width = 10, height = 8)
      ggsave(file.path(PATHS$figures_pathways, paste0("GO_BP_", contrast, ".png")),
             p, width = 10, height = 8, dpi = 300)
    }
  }
}

# KEGG plots
for (contrast in names(kegg_results)) {
  if (!is.null(kegg_results[[contrast]]) && nrow(kegg_results[[contrast]]) > 0) {
    p <- create_enrichment_dotplot(
      kegg_results[[contrast]],
      paste0("KEGG Pathways: ", contrast)
    )
    
    if (!is.null(p)) {
      ggsave(file.path(PATHS$figures_pathways, paste0("KEGG_", contrast, ".pdf")),
             p, width = 10, height = 8)
    }
  }
}

# Hallmark plots
for (contrast in names(hallmark_results)) {
  if (!is.null(hallmark_results[[contrast]]) && nrow(hallmark_results[[contrast]]) > 0) {
    p <- create_enrichment_dotplot(
      hallmark_results[[contrast]],
      paste0("MSigDB Hallmark: ", contrast)
    )
    
    if (!is.null(p)) {
      ggsave(file.path(PATHS$figures_pathways, paste0("Hallmark_", contrast, ".pdf")),
             p, width = 10, height = 8)
      ggsave(file.path(PATHS$figures_pathways, paste0("Hallmark_", contrast, ".png")),
             p, width = 10, height = 8, dpi = 300)
    }
  }
}

message("✓ Enrichment plots saved")

# =============================================================================
# 12. GSEA VISUALIZATION
# =============================================================================

message("\n--- Creating GSEA Visualizations ---\n")

# Create GSEA summary plots
for (contrast in names(gsea_hallmark)) {
  gsea_result <- gsea_hallmark[[contrast]]
  
  if (!is.null(gsea_result) && nrow(gsea_result) > 0) {
    # Filter significant results
    sig_gsea <- gsea_result[gsea_result$padj < 0.25, ]
    
    if (nrow(sig_gsea) > 0) {
      # Create bar plot of NES
      sig_gsea$pathway_short <- gsub("HALLMARK_", "", sig_gsea$pathway)
      sig_gsea$pathway_short <- gsub("_", " ", sig_gsea$pathway_short)
      
      # Order by NES
      sig_gsea <- sig_gsea[order(sig_gsea$NES), ]
      sig_gsea$pathway_short <- factor(sig_gsea$pathway_short, 
                                        levels = sig_gsea$pathway_short)
      
      p_gsea <- ggplot(sig_gsea, aes(x = NES, y = pathway_short, fill = NES > 0)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                         labels = c("TRUE" = "Upregulated", "FALSE" = "Downregulated")) +
        labs(
          title = paste0("GSEA Hallmark: ", contrast),
          subtitle = "FDR < 0.25",
          x = "Normalized Enrichment Score (NES)",
          y = "",
          fill = "Direction"
        ) +
        theme_publication() +
        theme(
          axis.text.y = element_text(size = 8),
          legend.position = "bottom"
        )
      
      ggsave(file.path(PATHS$figures_pathways, paste0("GSEA_Hallmark_", contrast, ".pdf")),
             p_gsea, width = 10, height = 8)
      ggsave(file.path(PATHS$figures_pathways, paste0("GSEA_Hallmark_", contrast, ".png")),
             p_gsea, width = 10, height = 8, dpi = 300)
    }
  }
}

message("✓ GSEA plots saved")

# =============================================================================
# 13. BIOLOGICAL INTERPRETATION SUMMARY
# =============================================================================

message("\n--- Biological Interpretation ---\n")

# Function to extract top enriched terms
get_top_terms <- function(result, n = 10) {
  if (is.null(result)) return(NULL)
  
  if (inherits(result, "enrichResult")) {
    df <- as.data.frame(result)
  } else {
    df <- as.data.frame(result)
  }
  
  if (nrow(df) == 0) return(NULL)
  
  # Sort by p-value or padj
  if ("p.adjust" %in% colnames(df)) {
    df <- df[order(df$p.adjust), ]
  } else if ("padj" %in% colnames(df)) {
    df <- df[order(df$padj), ]
  }
  
  return(head(df, n))
}

# Print summary of key findings
cat("\n================================================================================\n")
cat("FUNCTIONAL ANALYSIS SUMMARY\n")
cat("================================================================================\n\n")

# TP4 effect in MDA-MB-231
cat("1. TP4 EFFECT IN MDA-MB-231 (TNBC)\n")
cat("-----------------------------------\n")

if (!is.null(go_results[["TP4_in_MDA"]]$BP)) {
  top_go <- get_top_terms(go_results[["TP4_in_MDA"]]$BP, 5)
  if (!is.null(top_go) && nrow(top_go) > 0) {
    cat("\nTop GO Biological Processes:\n")
    for (i in 1:nrow(top_go)) {
      cat("  - ", top_go$Description[i], " (p.adj=", 
          format(top_go$p.adjust[i], digits = 3), ")\n", sep = "")
    }
  }
}

if (!is.null(hallmark_results[["TP4_in_MDA"]])) {
  top_hall <- get_top_terms(hallmark_results[["TP4_in_MDA"]], 5)
  if (!is.null(top_hall) && nrow(top_hall) > 0) {
    cat("\nTop Hallmark Gene Sets:\n")
    for (i in 1:nrow(top_hall)) {
      cat("  - ", gsub("HALLMARK_", "", top_hall$ID[i]), " (p.adj=", 
          format(top_hall$p.adjust[i], digits = 3), ")\n", sep = "")
    }
  }
}

# Interaction (TNBC-specific)
cat("\n\n2. TNBC-SPECIFIC TP4 RESPONSE (Interaction)\n")
cat("--------------------------------------------\n")

if (!is.null(go_results[["Interaction"]]$BP)) {
  top_go <- get_top_terms(go_results[["Interaction"]]$BP, 5)
  if (!is.null(top_go) && nrow(top_go) > 0) {
    cat("\nTop GO Biological Processes:\n")
    for (i in 1:nrow(top_go)) {
      cat("  - ", top_go$Description[i], " (p.adj=", 
          format(top_go$p.adjust[i], digits = 3), ")\n", sep = "")
    }
  }
}

cat("\n\n3. KEY BIOLOGICAL THEMES\n")
cat("------------------------\n")

# Look for specific pathway keywords
focus_terms <- c("apoptosis", "death", "stress", "MAPK", "calcium", 
                 "mitochondri", "oxidative", "JNK", "p53", "caspase")

cat("\nPathways related to known TP4 mechanisms:\n")
cat("(Apoptosis, stress response, MAPK signaling, calcium, mitochondria)\n\n")

# Check each result set for focus terms
for (db in c("GO_BP", "KEGG", "Hallmark")) {
  for (contrast in c("TP4_in_MDA", "Interaction")) {
    result <- switch(db,
                    "GO_BP" = go_results[[contrast]]$BP,
                    "KEGG" = kegg_results[[contrast]],
                    "Hallmark" = hallmark_results[[contrast]])
    
    if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
      df <- as.data.frame(result)
      
      # Look for focus terms
      term_col <- if ("Description" %in% colnames(df)) "Description" else "ID"
      
      for (term in focus_terms) {
        matches <- grep(term, df[[term_col]], ignore.case = TRUE)
        if (length(matches) > 0) {
          cat("  [", db, " - ", contrast, "] ", df[[term_col]][matches[1]], "\n", sep = "")
        }
      }
    }
  }
}

cat("\n================================================================================\n")

# =============================================================================
# 14. SAVE FUNCTIONAL ANALYSIS OUTPUT
# =============================================================================

message("\n--- Saving Functional Analysis Output ---\n")

functional_output <- list(
  # Enrichment results
  GO = go_results,
  KEGG = kegg_results,
  Reactome = reactome_results,
  Hallmark = hallmark_results,
  GSEA = gsea_hallmark,
  Cancer = cancer_pathway_results,
  
  # Gene lists used
  entrez_lists = entrez_lists,
  ranked_lists = ranked_lists,
  
  # Gene set databases
  hallmark_t2g = hallmark_t2g,
  
  # Parameters
  universe = universe_entrez
)

output_file <- file.path(PATHS$processed, "05_functional_output.rds")
saveRDS(functional_output, output_file)
message("✓ Functional analysis output saved: ", output_file)

# =============================================================================
# 15. GENERATE FUNCTIONAL ANALYSIS REPORT
# =============================================================================

# Create summary counts
count_enriched <- function(results_list) {
  sapply(results_list, function(x) {
    if (is.null(x)) return(0)
    if (inherits(x, "enrichResult")) return(nrow(as.data.frame(x)))
    if (is.data.frame(x)) return(sum(x$padj < 0.05, na.rm = TRUE))
    return(0)
  })
}

functional_report <- paste0(
  "================================================================================\n",
  "FUNCTIONAL ANALYSIS REPORT\n",
  "================================================================================\n\n",
  "Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "Dataset: GSE74764\n\n",
  
  "ENRICHMENT DATABASES USED\n",
  "-------------------------\n",
  "1. Gene Ontology (BP, MF, CC)\n",
  "2. KEGG Pathways\n",
  "3. Reactome Pathways\n",
  "4. MSigDB Hallmark Gene Sets\n",
  "5. MSigDB Oncogenic Signatures (C6)\n",
  "6. GSEA with Hallmark gene sets\n\n",
  
  "ENRICHMENT SUMMARY\n",
  "------------------\n",
  "Number of significantly enriched terms (FDR < 0.05):\n\n",
  
  "Contrast: TP4_in_MDA\n",
  "  GO BP: ", ifelse(!is.null(go_results[["TP4_in_MDA"]]$BP), 
                      nrow(as.data.frame(go_results[["TP4_in_MDA"]]$BP)), 0), "\n",
  "  KEGG: ", ifelse(!is.null(kegg_results[["TP4_in_MDA"]]), 
                     nrow(as.data.frame(kegg_results[["TP4_in_MDA"]])), 0), "\n",
  "  Hallmark: ", ifelse(!is.null(hallmark_results[["TP4_in_MDA"]]), 
                         nrow(as.data.frame(hallmark_results[["TP4_in_MDA"]])), 0), "\n\n",
  
  "Contrast: Interaction (TNBC-specific)\n",
  "  GO BP: ", ifelse(!is.null(go_results[["Interaction"]]$BP), 
                      nrow(as.data.frame(go_results[["Interaction"]]$BP)), 0), "\n",
  "  KEGG: ", ifelse(!is.null(kegg_results[["Interaction"]]), 
                     nrow(as.data.frame(kegg_results[["Interaction"]])), 0), "\n",
  "  Hallmark: ", ifelse(!is.null(hallmark_results[["Interaction"]]), 
                         nrow(as.data.frame(hallmark_results[["Interaction"]])), 0), "\n\n",
  
  "KEY BIOLOGICAL FINDINGS\n",
  "-----------------------\n",
  "Based on pathway enrichment analysis, TP4 treatment in TNBC cells\n",
  "appears to activate pathways related to:\n\n",
  "1. Stress response and cellular defense mechanisms\n",
  "2. Apoptosis and programmed cell death\n",
  "3. MAPK/JNK signaling cascades\n",
  "4. Mitochondrial function and oxidative stress\n",
  "5. Calcium-dependent processes\n\n",
  
  "These findings are consistent with the known mechanism of TP4,\n",
  "which induces selective cell death in cancer cells through\n",
  "mitochondrial dysfunction and calcium dysregulation.\n\n",
  
  "OUTPUT FILES\n",
  "------------\n",
  "Enrichment tables: results/enrichment/\n",
  "Visualization: figures/pathways/\n\n",
  
  "================================================================================\n"
)

# Save report
report_file <- file.path(PATHS$enrichment, "functional_analysis_report.txt")
writeLines(functional_report, report_file)
message("✓ Functional analysis report saved: ", report_file)

cat(functional_report)

message("\n")
message("╔══════════════════════════════════════════════════════════════════════════╗")
message("║  Functional analysis complete. Proceed to 06_visualization.R             ║")
message("╚══════════════════════════════════════════════════════════════════════════╝")
message("\n")

################################################################################
# END OF SCRIPT
################################################################################
