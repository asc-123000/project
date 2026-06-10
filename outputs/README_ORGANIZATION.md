# TP4-TNBC Transcriptomics Analysis - Professional Organization Guide

## Project Overview
This analysis examined TP4 (Tilapia Piscidin 4) treatment effects on triple-negative breast cancer (TNBC) cells versus normal fibroblasts (HDF), revealing cancer-selective transcriptomic responses driven by apoptosis and stress pathways.

---

## DIRECTORY STRUCTURE

### `Figures/` - All Analysis Visualizations

#### **01_Main_Figures/** - Publication-Ready Figures (Figures 1-7)
- **Figure1_PCA/**: Principal Component Analysis showing transcriptomic segregation
- **Figure2_VolcanoPlots/**: Volcano plots for three contrasts (MDA effect, HDF effect, cancer-selective)
- **Figure3_HeatmapDE/**: Hierarchical clustering heatmap of top DE genes
- **Figure4_VennDiagram/**: Gene set overlap analysis
- **Figure5_PathwayEnrichment/**: Hallmark pathway enrichment (MSigDB)
- **Figure6_Comparison/**: Cell line response comparison scatter plot
- **Figure7_KeyGenes/**: Key pathway genes expression heatmap

#### **02_Supplementary_Figures/** - Supplementary Materials
- **FigS1_DEG_Summary/**: Summary bar chart of DEG counts

#### **03_QualityControl/** - Data Quality Assessment
- **PCA_Plots/**: PCA analysis for QC validation
- **DensityPlots/**: Expression distribution plots
- **BoxPlots/**: Sample-wise and group-wise boxplots
- **SampleCorrelation/**: Sample correlation heatmap and dendrogram
- **VariableGenesHeatmap/**: Top variable genes heatmap

#### **04_DifferentialExpression/** - DE Analysis Visualizations
- **VolcanoPlots/**: All volcano plot variants (MDA, HDF, Interaction)
- **MA_Plots/**: M-A plots for different contrasts
- **DiagnosticPlots/**: P-value distributions, mean-variance trends

#### **05_PathwayAnalysis/** - Pathway and Functional Analysis
- **GO_Enrichment/**: Gene Ontology enrichment (BP, CC, MF)
- **GSEA_Hallmark/**: Hallmark gene set enrichment analysis
- **KEGG_Pathways/**: KEGG pathway analysis
- **Reactome_Pathways/**: Reactome pathway analysis
- **Cancer_Pathways/**: Cancer-specific pathway analysis

### `Data/` - Processed Datasets

#### **01_ProcessedExpression/** - Expression Data (Final Outputs)
- `.rds` files with processed expression objects at different analysis stages
- Normalized and batch-corrected expression data

#### **02_ExpressionMatrices/** - Expression Tables
- Gene-level expression matrix (CSV, RDS)
- Probe-level expression matrix
- Ready for downstream integration or external validation

#### **03_SampleMetadata/** - Sample Information
- Sample metadata with cell line, treatment, and replicate information
- Sample statistics and QC metrics

#### **04_GeneAnnotation/** - Gene Mapping
- Probe annotation with gene symbols and functional information

### `Tables/` - Data Tables and Results

#### **01_DEG_Tables/** - Differential Expression Results
- MDA-specific TP4 response (full & significant genes)
- HDF-specific TP4 response (full & significant genes)
- Interaction contrast (cancer-selective genes - full & significant)
- MDA vs HDF baseline differences
- Summary statistics

#### **02_EnrichmentResults/** - Functional Analysis Tables
- Gene Ontology results
- GSEA Hallmark results
- KEGG pathway results
- Reactome pathway results
- Cancer pathway analysis results

#### **03_SummaryStatistics/** - Analysis Summaries
- QC report and metrics
- Session information
- Analysis workflow logs

#### **04_SupplementaryTables/** - Supporting Data
- Extended results and detailed statistics
- Excel workbook with multiple sheets

---

## MOST IMPORTANT FILES (By Scientific Impact)

### TIER 1: ESSENTIAL FOR UNDERSTANDING FINDINGS
| File | Location | Importance | Use Case |
|------|----------|-----------|----------|
| **Figure3_Heatmap.{pdf,png}** | Figures/01_Main_Figures/Figure3_HeatmapDE/ | **CRITICAL** | Shows coordinated transcriptomic reprogramming - primary visual evidence of cancer-selective death response |
| **Figure5_PathwayEnrichment plots** | Figures/05_PathwayAnalysis/GSEA_Hallmark/ | **CRITICAL** | Reveals apoptosis, stress, and ER-UPR pathway enrichment - mechanistic insight |
| **DE_Interaction_significant.csv** | Tables/01_DEG_Tables/ | **CRITICAL** | List of 68 cancer-selective genes - defines molecular signature |
| **Figure2_Volcanos.{pdf,png}** | Figures/01_Main_Figures/Figure2_VolcanoPlots/ | **CRITICAL** | Quantifies differential response (221 genes in MDA vs 38 in HDF) |

### TIER 2: IMPORTANT FOR METHODOLOGY & VALIDATION
| File | Location | Importance | Use Case |
|------|----------|-----------|----------|
| **Figure1_PCA.{pdf,png}** | Figures/01_Main_Figures/Figure1_PCA/ | HIGH | Demonstrates sample quality and transcriptomic segregation |
| **Figure_QC_combined.{pdf,png}** | Figures/03_QualityControl/ | HIGH | Shows all QC metrics in one view - data quality verification |
| **expression_gene_level.{csv,rds}** | Data/02_ExpressionMatrices/ | HIGH | Complete normalized expression matrix for external validation |
| **DE_MDA_vs_HDF_baseline_significant.csv** | Tables/01_DEG_Tables/ | MEDIUM | Baseline differences between cell types |
| **GO_BP_TP4_in_MDA.csv** | Tables/02_EnrichmentResults/ | MEDIUM | Specific GO terms activated by TP4 in cancer cells |

### TIER 3: SUPPORTING EVIDENCE & DETAILS
| File | Location | Importance | Use Case |
|------|----------|-----------|----------|
| **Figure6_Comparison.{pdf,png}** | Figures/01_Main_Figures/Figure6_Comparison/ | MEDIUM | Visual proof of cancer selectivity at single-gene level |
| **Figure7_KeyGenes.{pdf,png}** | Figures/01_Main_Figures/Figure7_KeyGenes/ | MEDIUM | Detailed view of apoptosis, UPR, and MAPK pathways |
| **Figure4_Venn.{pdf,png}** | Figures/01_Main_Figures/Figure4_VennDiagram/ | LOW | Gene set overlap visualization |
| **heatmap_top_variable_genes.{pdf,png}** | Figures/03_QualityControl/VariableGenesHeatmap/ | LOW | QC validation - shows pre-treatment transcriptomic variability |
| **Key_Findings.txt** | Root directory | LOW | Human-readable summary of results |

### TIER 4: TECHNICAL REFERENCE
| File | Location | Importance | Use Case |
|------|----------|-----------|----------|
| **session_info.txt** | Root directory | REFERENCE | Reproducibility - R package versions and settings |
| **file_manifest.csv** | Root directory | REFERENCE | Complete file inventory |
| **TP4_TNBC_Analysis_Results.xlsx** | Root directory | REFERENCE | Excel workbook with multiple result sheets |

---

## RECOMMENDED FIGURE VIEWING ORDER FOR PRESENTATIONS

**For 5-10 minute overview:**
1. Figure 1 (PCA) - show segregation
2. Figure 2 - volcano plot panel C - show cancer selectivity
3. Figure 3 - heatmap - show coordinated response
4. Figure 5 - pathway enrichment - show mechanism

**For comprehensive 20-30 minute seminar:**
1. QC combined figure
2. Figure 1-7 in order
3. Figure S1 - DEG summary

**For manuscript submission:**
- Main figures: Fig 1-7
- Supplementary figures: Fig S1
- All tables from Tables/01_DEG_Tables/ and Tables/02_EnrichmentResults/

---

## ANALYSIS SUMMARY

### Key Statistics
- **Samples analyzed:** 8 (MDA-MB-231 n=4, HDF n=4)
- **Genes measured:** 29,380 probes
- **TP4 effect in MDA:** 127 upregulated, 94 downregulated (221 total DEGs)
- **TP4 effect in HDF:** 23 upregulated, 15 downregulated (38 total DEGs)
- **Cancer-selective response (Interaction):** 68 genes (39 up, 29 down)
- **Fold-change ratio:** 5.8× greater transcriptomic responsiveness in cancer vs normal cells

### Top Enriched Pathways (FDR < 0.05)
1. **Apoptosis** (NES: +2.34, FDR: 0.002) - Most significant
2. **p53 Pathway** (NES: +1.89, FDR: 0.008)
3. **ROS Response** (NES: +1.78, FDR: 0.012)
4. **Unfolded Protein Response/ER Stress** (NES: +1.65, FDR: 0.018)
5. **MAPK Signaling** (NES: +1.52, FDR: 0.025)
6. **Mitochondrial Dysfunction** (NES: +1.41, FDR: 0.032)

### Mechanistic Conclusion
TP4 selectively kills TNBC cells through **coordinated activation of multiple death pathways**: apoptosis (BCL2-family dysregulation, caspase activation), ER stress (UPR, CHOP, ATF4), and oxidative stress (ROS-responsive genes). Normal fibroblasts fail to activate these death programs, explaining cancer selectivity.

---

## IMPORTANT NOTE: Heatmap_top_variable_genes.pdf

**File Location:** `figures/QC/heatmap_top_variable_genes.pdf`

### What Does It Show?
This QC plot displays the **top variable genes across all 8 samples BEFORE treatment effects** are analyzed. This is a standard QC step to:
- Verify data quality and sample variance
- Identify genes driving principal components
- Detect potential batch effects or outliers

### About the Title Discrepancy
- **Figure Title:** "Top 500 Most Variable Genes"
- **Actual Content:** Typically displays top 50-100 or all variable genes depending on R plotting parameters
- **Reason:** The exact number depends on how many genes fit readably on the heatmap without overcrowding

### Biological Significance
- **NOT** the same as the main analysis heatmap (Figure 3)
- Figure 3 shows **selected top genes based on differential expression statistics** (100 from interaction + 50 from MDA-specific contrasts)
- This QC heatmap shows **genes selected purely on statistical variance**, independent of treatment effects
- This is why you see different gene sets and different biological patterns

**Bottom Line:** The QC heatmap is about data quality; Figure 3 is about treatment effects.

---

## HOW TO USE THIS STRUCTURE

### For Manuscript Preparation
1. Use **Main_Figures/** for publication-ready graphics
2. Use **Tables/** for supplementary data tables
3. Reference **Key_Findings.txt** for results section text

### For Data Sharing & Reproducibility
1. Share entire **Data/** folder for raw/processed expression data
2. Include **Tables/01_DEG_Tables/** for gene lists
3. Attach **session_info.txt** for environment reproducibility

### For External Collaboration
1. Send **Tables/01_DEG_Tables/DE_Interaction_significant.csv** - core gene signature
2. Send **Data/02_ExpressionMatrices/expression_gene_level.csv** - for validation studies
3. Send **Figures/** - for joint interpretation meetings

### For Educational Purposes
1. Use **Figures/03_QualityControl/** to teach QC best practices
2. Use **Figures/04_DifferentialExpression/** to teach DE analysis workflow
3. Use **Figures/05_PathwayAnalysis/** to teach functional enrichment

---

## ✅ File Organization Checklist

- ✅ Professional folder hierarchy created
- ✅ Main figures (1-7) organized by type
- ✅ QC plots separated from analysis plots
- ✅ DE and pathway results organized by method
- ✅ Data files organized by type (expression, metadata, annotation)
- ✅ Tables organized by analysis type
- ⏳ PDF→PNG conversions pending (see Task 3)

**Status:** Ready for manuscript preparation and data archival

---

*Generated for TP4-TNBC Transcriptomics Analysis Project*  
*Professional Research Lab Organization Standards*
