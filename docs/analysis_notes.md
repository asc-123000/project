# Analysis Notes: TP4 TNBC Transcriptomics

## Technical Decisions and Scientific Justifications

### 1. Using Pre-Processed Data vs Raw Data

**Decision**: Use GEO series matrix (pre-processed) rather than raw Agilent files.

**Justification**:
- With n=2 per group, we have minimal tolerance for additional technical variability
- Raw data preprocessing requires choices (background correction, normalization) that could introduce variance
- Original authors applied appropriate Agilent-specific processing
- QC verification confirms adequate normalization

**Trade-off**: Less control over preprocessing, but more stable results given sample size.

---

### 2. Statistical Thresholds

#### FDR < 0.05
- Standard discovery threshold
- Appropriate for microarray experiments
- Balances false positives and false negatives

#### |log2FC| > 0.585 (1.5-fold)
- More lenient than typical 2-fold cutoff
- Justified by:
  - Small sample size limiting power
  - 6-hour treatment may show subtle effects
  - Biologically meaningful changes can be < 2-fold

---

### 3. Why Cell-Means Model?

**Alternative**: Could use factorial model with main effects + interaction:
```r
design <- model.matrix(~ cell_line * treatment, data = sample_metadata)
```

**Our choice**: Cell-means model (group model):
```r
design <- model.matrix(~ 0 + group, data = sample_metadata)
```

**Rationale**:
- More intuitive contrast specification
- Direct estimation of group means
- Clearer biological interpretation
- Same statistical results, different parameterization

---

### 4. Interaction Contrast Interpretation

The interaction contrast is:
```
(MDA_TP4 - MDA_Mock) - (HDF_TP4 - HDF_Mock)
```

This identifies genes where the TP4 effect differs between cell types.

**Example interpretations**:

| Gene | MDA response | HDF response | Interaction | Interpretation |
|------|--------------|--------------|-------------|----------------|
| GeneA | +3 | +1 | +2 | More induced in cancer |
| GeneB | -2 | 0 | -2 | Suppressed only in cancer |
| GeneC | +2 | +2 | 0 | Common response (not specific) |
| GeneD | +1 | +3 | -2 | Less induced in cancer |

---

### 5. Empirical Bayes Rationale

With n=2 per group, gene-specific variance estimates are unreliable.

limma's eBayes:
- Borrows information across ~20,000 genes
- Stabilizes variance estimates
- Shrinks extreme variances toward a common prior
- Dramatically improves power for small n

**Evidence**: Prior degrees of freedom (df.prior) indicates degree of shrinkage.
- df.prior > 3 suggests substantial information borrowing
- df.prior < 1 suggests limited shrinkage (poor gene-wise variance estimation)

---

### 6. Probe-to-Gene Mapping Strategy

**Problem**: Multiple probes can map to the same gene.

**Strategy**: Keep the probe with highest mean expression.

**Rationale**:
- Highest-expressed probe likely represents the most abundant/relevant transcript
- Avoids artificial variance inflation from averaging
- Common approach in microarray analysis

**Alternative considered**: Average across probes (rejected due to potential noise from low-expression probes).

---

### 7. Pathway Analysis Databases

| Database | Purpose | Why Included |
|----------|---------|--------------|
| GO (BP, MF, CC) | General biological processes | Comprehensive coverage |
| KEGG | Metabolic and signaling pathways | Well-curated, pathway-centric |
| Reactome | Detailed molecular reactions | Mechanistic detail |
| Hallmark | Curated, coherent gene sets | Reduces redundancy |
| GSEA | Ranked enrichment | Uses all genes, not just significant |

**Excluded**: TCGA expression signatures (not comparable to microarray data).

---

### 8. Focus Pathways

Based on TP4 mechanism (Ting et al., 2016):

1. **Apoptosis/Cell Death**
   - TP4 induces selective cell death
   - Look for: BCL2 family, caspases, death receptors

2. **ER Stress/UPR**
   - TP4 causes calcium release from ER
   - Look for: ATF4, DDIT3 (CHOP), XBP1, HSPA5 (BiP)

3. **MAPK/JNK/AP-1**
   - Original study highlighted FOSB
   - Look for: FOS, JUN family, MAPK pathway

4. **Calcium Signaling**
   - TP4 disrupts calcium homeostasis
   - Look for: calcium-dependent enzymes

5. **Mitochondrial Dysfunction**
   - TP4 affects mitochondrial membrane potential
   - Look for: OXPHOS, ROS response

---

### 9. Figure Standards

All figures designed for journal submission:
- 300 DPI minimum
- PDF (vector) for line art
- PNG for complex heatmaps
- Consistent color scheme
- Clear, minimal design
- Proper axis labels and legends

---

### 10. Known Limitations

1. **Sample Size**: n=2 severely limits power
2. **Single Time Point**: 6h may miss early/late responses
3. **Single Cell Line**: MDA-MB-231 may not represent all TNBC
4. **In Vitro**: Lacks tumor microenvironment
5. **Pre-processed Data**: Limited QC control

---

## Troubleshooting

### Package Installation Issues

```r
# If BiocManager fails:
options(repos = BiocManager::repositories())
install.packages("BiocManager")
BiocManager::install(version = "3.17")  # or current version

# If specific package fails:
BiocManager::install("package_name", force = TRUE)
```

### GEO Download Issues

```r
# If GEOquery times out:
options(timeout = 300)  # Increase timeout

# Or download manually from GEO website
```

### Memory Issues

```r
# For large operations:
gc()  # Force garbage collection

# Reduce plot resolution if needed
ANALYSIS_PARAMS$plot_dpi <- 150
```

---

*Notes compiled during analysis development*
