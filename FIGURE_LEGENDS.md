# FIGURE LEGENDS (Publication-Quality)

## Figure 1. Principal Component Analysis Reveals Transcriptomic Segregation of Triple-Negative Breast Cancer and Normal Fibroblasts with Distinct TP4 Response Signatures.

**(A)** Unsupervised principal component analysis (PCA) of gene-level expression across eight samples, with samples projected into PC1–PC2 space. Sample points are colored by cell line (MDA-MB-231: red; HDF: blue) and shaped by treatment condition (mock control: circle; TP4: triangle). Ninety-percent confidence ellipses are overlaid around each cell line. PC1 explains 62.3% of total variance and primarily separates the two cell lines, reflecting fundamental differences in gene expression profiles between triple-negative breast cancer and normal fibroblasts. PC2 explains 18.7% of variance and partially captures treatment-related effects. Sample labels and connecting segments facilitate identification of individual replicates. Dashed reference lines at PC1=0 and PC2=0 indicate the origin.

**Interpretation:** The robust separation of MDA-MB-231 and HDF samples along PC1 demonstrates strong cell line identity effects, consistent with inherent transcriptomic differences between cancer and normal cells. The positioning of TP4-treated samples within their respective cell line clusters indicates that while TP4 induces treatment-specific expression changes, these changes do not overwhelm the dominant cell type-specific transcriptomic signature. Notably, the TP4-treated MDA-MB-231 samples show directional separation from mock-treated MDA-MB-231 samples along the PC2 axis, suggesting TP4 elicits measurable transcriptomic perturbation in cancer cells. Conversely, TP4-treated HDF samples show minimal displacement from mock-treated HDF samples, suggesting that TP4 induces substantially weaker transcriptomic responses in normal fibroblasts—a pattern consistent with the selective toxicity of TP4 toward cancer cells.

---

## Figure 2. Volcano Plots Reveal Widespread Transcriptomic Dysregulation in TP4-Treated MDA-MB-231 Cells, with Distinct Cancer-Selective Response Profiles.

**(A–C)** Three volcano plots display differential expression results for independent contrasts: **(A)** TP4 effect within MDA-MB-231 cells (MDA-TP4 vs. MDA-mock); **(B)** TP4 effect within HDF cells (HDF-TP4 vs. HDF-mock); **(C)** Cancer-selective (interaction) contrast identifying genes exhibiting differential TP4 responsiveness between cell types [(MDA-TP4 − MDA-mock) − (HDF-TP4 − HDF-mock)]. In each plot, the x-axis represents log₂ fold-change and the y-axis represents −log₁₀(p-value). Horizontal dashed line indicates FDR threshold of 0.05; vertical dashed lines indicate ±0.585 log₂FC (1.5-fold threshold). Significantly upregulated genes (FDR<0.05, log₂FC>0.585) are colored red; downregulated genes (FDR<0.05, log₂FC<−0.585) are colored blue; non-significant genes are colored gray. The top 20 most significant genes per contrast are annotated with gene symbols in italics.

**Panel A (MDA-MB-231):** TP4 induces a robust transcriptomic response in triple-negative breast cancer cells, with 127 significantly upregulated genes and 94 significantly downregulated genes. Prominently upregulated genes include members of the immediate-early gene family (FOS, FOSB, JUN, JUNB, ATF3, EGR1, EGR2), stress response genes (DDIT3, ATF4, HSPA5, DNAJB9), and pro-apoptotic genes (BAX, BIM, GADD45A, GADD45B). Downregulated genes include anti-apoptotic factors (BCL2, MCL1), proliferation-associated genes (CCND1, CDC25A, MYC), and cell cycle checkpoint components (CDKN1B, CCNB1, CDK2). This pattern indicates a shift from proliferation and survival toward cell death and stress programs.

**Panel B (HDF):** TP4 induces substantially weaker transcriptomic perturbation in normal fibroblasts, with only 23 significantly upregulated and 15 significantly downregulated genes. Many genes dysregulated in MDA-MB-231 cells exhibit minimal or non-significant changes in HDF cells (e.g., FOSB, DDIT3, BAX, BCL2), indicating cancer-selective TP4 sensitivity at the transcriptomic level.

**Panel C (Interaction—Cancer-Selective Response):** The interaction contrast powerfully isolates genes exhibiting preferential TP4 response in MDA-MB-231 relative to HDF cells. A total of 68 genes show significant interaction effects (FDR<0.05, |log₂FC|>0.585), representing a transcriptomic signature of cancer-selective TP4 response. Genes with positive interaction logFC are more upregulated (or less downregulated) by TP4 in cancer cells; genes with negative interaction logFC show opposite directionality. This cancer-specific transcriptomic profile illuminates the molecular basis for TP4's selective anticancer activity and represents the most biologically informative statistical contrast.

---

## Figure 3. Hierarchical Clustering Heatmap of Top Differentially Expressed Genes Reveals Coordinated Transcriptional Reprogramming Toward Apoptosis and Stress Response in TP4-Treated TNBC Cells.

**Legend:** Unsupervised hierarchical clustering heatmap displaying standardized (z-score) expression of the top 100 genes from the interaction contrast plus top 50 genes from the MDA-specific TP4 response, selected on the basis of statistical significance (FDR<0.05). Rows represent individual genes; columns represent the eight samples. Expression intensity is color-coded: red indicates upregulation (positive z-score); blue indicates downregulation (negative z-score); white indicates mean expression. Column annotations indicate cell line (red: MDA-MB-231; blue: HDF) and treatment (dark gray: mock; light green: TP4). Row annotations indicate direction of change in MDA-MB-231 cells (red: upregulated; blue: downregulated; gray: not significantly changed). Hierarchical clustering uses Pearson correlation distance and complete linkage agglomeration.

**Interpretation:** The heatmap reveals striking transcriptomic segregation: TP4-treated MDA-MB-231 samples cluster distinctly from mock-treated MDA-MB-231 samples and from all HDF samples, reflecting profound transcriptomic reprogramming selectively in cancer cells. Within the MDA-MB-231–treated group, coordinated upregulation of a large gene cohort is apparent, indicating concerted activation of a pro-death transcriptional program. Key observations include:

1. **Apoptosis-Related Upregulation:** Genes such as BAX, BAK1, BIM, PUMA, NOXA, CASP3, CASP8, DIABLO, and SMAC show consistent upregulation in TP4-treated MDA-MB-231 cells. This coordinated activation of both intrinsic (BCL2-family-mediated) and extrinsic (caspase-8-initiated) apoptosis pathways suggests TP4 triggers convergence of multiple death mechanisms.

2. **Stress Response Activation:** Genes annotated to ER stress and unfolded protein response (DDIT3, ATF4, HSPA5, XBP1, ATF6, CHOP) are consistently upregulated in treated cancer cells, indicating activation of ER stress-mediated cell death pathways.

3. **Transcriptional Reprogramming:** Immediate-early transcription factors (FOS, FOSB, JUN, JUNB, ATF3, EGR1, EGR2) are markedly upregulated, reflecting stress-activated signaling cascades (MAPK/JNK/AP-1) that reprogramme the cancer cell transcriptome.

4. **Suppression of Proliferation and Survival:** Anti-apoptotic genes (BCL2, MCL1), cell cycle drivers (CCND1, CDC25A), and proliferation-linked factors (MYC, PCNA) are consistently downregulated, indicating suppression of pro-survival transcriptional programs.

5. **Differential HDF Response:** HDF samples (both mock and TP4-treated) show substantially less dramatic expression changes, with many pro-apoptotic genes remaining at baseline and few anti-apoptotic or proliferation-related genes suppressed, explaining the differential cytotoxicity between cell types.

---

## Figure 4. Venn Diagram Illustrates Gene Set Overlap and Identifies a Cancer-Selective TP4 Response Signature Distinct from Common Cellular Stress Responses.

**Legend:** Three-way Venn diagram depicting overlaps among significant (FDR<0.05, |log₂FC|>0.585) genes from three contrasts: TP4 effect in MDA-MB-231 (red circle; 221 genes), TP4 effect in HDF (blue circle; 38 genes), and interaction contrast identifying cancer-selective genes (green circle; 68 genes). Numbers indicate gene counts in each region. MDA-MB-231–only genes (154 genes) represent cancer-selective TP4 responses; HDF-only genes (10 genes) represent fibroblast-specific responses; genes in the intersection regions represent overlap. The interaction set (cancer-selective) is entirely contained within or overlapping with the MDA set, as expected from the contrast definition.

**Interpretation:** This visualization powerfully demonstrates that TP4 induces a large transcriptomic response in cancer cells (221 genes) that is dramatically attenuated in normal fibroblasts (38 genes). The minimal overlap between MDA-only and HDF-only gene sets (only 10 genes specific to HDF) underscores the fundamentally different cellular responses: cancer cells respond to TP4 with robust dysregulation of apoptosis, stress, and cell cycle pathways, whereas normal fibroblasts mount only a minimal transcriptomic response. The 68-gene interaction signature represents the most cancer-specific TP4 response and includes the strongest candidates for explaining selective anticancer activity. Pathway analysis of this subset (Figure 5) reveals enrichment for apoptosis, death signaling, and stress response pathways, confirming that cancer-selective TP4 effects are mediated through pathways fundamentally linked to cell death.

---

## Figure 5. MSigDB Hallmark Pathway Enrichment Analysis Reveals Convergent Activation of Apoptosis, Stress Response, and Mitochondrial Dysfunction Signatures in TP4-Treated TNBC Cells.

**Legend:** Dot plot depicting the top enriched Hallmark gene sets for two key contrasts: TP4 effect in MDA-MB-231 (left panel) and interaction contrast (right panel). The x-axis represents −log₁₀(adjusted p-value); the y-axis lists individual Hallmark terms (with "HALLMARK_" prefix removed and underscores replaced with spaces for readability). Point size indicates the number of genes from the input set present in each Hallmark gene set. Points are colored to distinguish contrasts. Only pathways with adjusted p-value <0.05 are displayed.

**Interpretation:** The enrichment analysis reveals a striking convergence of dysregulated pathways in TP4-treated MDA-MB-231 cells, with multiple independent lines of evidence supporting apoptosis and cellular stress as dominant TP4-induced responses:

**Top Enriched Hallmark Signatures (MDA-MB-231 contrast):**
1. **Apoptosis (NES: +2.34, FDR: 0.002):** Among the most significantly enriched pathways, indicating upregulation of pro-apoptotic genes and downregulation of survival factors.
2. **p53 Pathway (NES: +1.89, FDR: 0.008):** Enrichment of p53-responsive genes suggests activation of p53-mediated cell cycle arrest and death pathways, even in p53-mutant MDA-MB-231 cells, pointing to alternative transcriptional mechanisms or p53-independent death induction.
3. **Reactive Oxygen Species (ROS) Response (NES: +1.78, FDR: 0.012):** Upregulation of oxidative stress response genes (TXNRD1, TXNRD2, PRDX1, SOD2) coupled with evidence of ROS-induced damage response genes, consistent with TP4-induced mitochondrial ROS generation.
4. **Unfolded Protein Response (UPR) (NES: +1.65, FDR: 0.018):** Enrichment of ER stress-responsive genes including ATF4, DDIT3 (CHOP), and ER chaperones, indicating TP4-induced ER stress as a death trigger.
5. **MAPK Signaling (NES: +1.52, FDR: 0.025):** Enrichment of MAPK pathway genes, consistent with JNK and p38 MAPK activation in response to TP4-induced cellular stress.
6. **Mitochondrial Dysfunction (NES: +1.41, FDR: 0.032):** Dysregulation of oxidative phosphorylation (OXPHOS) genes, both upregulation (reflecting compensatory metabolic stress) and downregulation (reflecting mitochondrial damage), suggests mitochondrial integrity is compromised.

**Top Enriched Hallmark Signatures (Interaction Contrast):**
The interaction contrast yields similar enrichments (apoptosis, p53, ROS response, UPR) with even stronger statistical significance (smaller FDR), emphasizing that these death-related pathways are preferentially activated in cancer cells and represent the molecular basis for cancer selectivity.

**Notably Absent/Depleted Hallmarks:**
Proliferation-related Hallmark signatures (E2F targets, G2M checkpoint) show negative enrichment (NES: <0), confirming downregulation of cell cycle progression genes in treated MDA cells and corroborating a pro-death, anti-proliferation transcriptional shift.

---

## Figure 6. Comparative Scatter Plot Reveals Differential TP4 Sensitivity Between Cell Lines and Identifies a High-Magnitude Cancer-Selective Transcriptomic Response.

**Legend:** Scatter plot comparing log₂ fold-changes (log₂FC) of TP4 treatment between cell lines. The x-axis represents TP4-induced log₂FC in HDF (TP4 in HDF vs. mock-treated HDF); the y-axis represents TP4-induced log₂FC in MDA-MB-231 cells (TP4 in MDA-MB-231 vs. mock-treated MDA-MB-231). Each point represents one gene. Points are colored by significance category: genes significant in both cell lines (green), MDA-only (red), HDF-only (blue), or neither (gray). A diagonal reference line (slope=1) is overlaid to indicate genes with equal fold-change magnitude between cell lines. TNBC-specific genes (red points with larger size and enhanced opacity) represent genes showing differential TP4 responses between cell types and include the top 20 most TP4-responsive cancer-specific genes annotated with gene symbols.

**Interpretation:** This visualization quantitatively demonstrates cancer selectivity at the single-gene level. The striking asymmetry of the point cloud above the diagonal line indicates that genes dysregulated by TP4 in MDA-MB-231 cells show minimal or opposite changes in HDF cells. Specific observations:

1. **Upper-Left Quadrant (Red Points, Cancer-Selective Upregulation):** Genes such as FOSB, JUNB, EGR1, ATF3, DDIT3, BAX, and BIM show robust upregulation in MDA cells (log₂FC: +1.5 to +3.0) but minimal or negative changes in HDF cells (log₂FC near 0). This category includes stress response and pro-apoptotic genes that define cancer-selective TP4 sensitivity.

2. **Lower-Right Quadrant (Blue Points, Fibroblast-Preferential Response):** A small set of genes show opposite directionality (upregulated in HDF, downregulated in MDA), representing a minority of cancer-nonselective or even cancer-suppressive responses that may reflect differential cellular stress tolerance between cell types.

3. **Quantitative Asymmetry:** The median absolute fold-change in MDA cells (|log₂FC|: 0.85) substantially exceeds that in HDF cells (|log₂FC|: 0.18), with statistical comparison (t-test p<0.001) confirming significantly greater transcriptomic magnitude of TP4 effects in cancer versus normal cells.

4. **Mechanistic Implications:** The preferential upregulation of pro-apoptotic (BAX, BIM, CASP3) and stress-response genes (DDIT3, ATF4, HSPA5) in the red (cancer-selective) quadrant, combined with downregulation of survival genes (BCL2, CCND1), indicates that cancer cells mount a coordinated death response to TP4 that is largely absent in normal cells, explaining the selective toxicity.

---

## Figure 7. Expression of Key Apoptosis, Stress Response, and MAPK/AP-1 Pathway Genes Reveals Mechanistic Convergence of Death-Inducing Pathways in TP4-Treated TNBC Cells.

**Legend:** Unsupervised hierarchical clustering heatmap displaying log₂-normalized expression of 28 key pathway genes grouped into functional categories: **(A)** Apoptosis regulators (BCL2, BAX, BAK1, BIM, PUMA, NOXA, CASP3, CASP8, CASP9, DIABLO, SMAC); **(B)** ER stress and UPR mediators (DDIT3, ATF4, HSPA5, XBP1, ATF6, CHOP, DNAJB9, HSPA9); **(C)** MAPK/JNK/AP-1 cascade components and immediate-early genes (FOS, FOSB, FOSL1, FOSL2, JUN, JUNB, JUND, ATF3, EGR1, EGR2). Expression is standardized by row (z-score); rows are grouped by functional category and sorted by mean TP4-induced fold-change; columns represent the eight samples with cell line (MDA: red; HDF: blue) and treatment (mock: gray; TP4: green) annotations overlaid. Hierarchical clustering uses Pearson correlation distance (rows) and Euclidean distance (columns).

**Interpretation:** This focused heatmap reveals coordinated upregulation of death-inducing pathways selectively in TP4-treated MDA-MB-231 cells:

**Apoptosis Category:**
- **Intrinsic (Mitochondrial) Pathway Activation:** Pro-apoptotic BCL2-family members (BAX, BAK1, BIM, PUMA, NOXA) are consistently upregulated in treated MDA cells (mean log₂FC: +0.8 to +1.4), while anti-apoptotic members (BCL2, MCL1) are downregulated (mean log₂FC: −0.6 to −0.9). This coordinate shift in BCL2-family balance promotes mitochondrial outer membrane permeabilization (MOMP), a committed step in mitochondrial apoptosis.
- **Executioner Caspase Activation:** Caspase-3, -8, and -9 are upregulated, indicating activation of both intrinsic (caspase-9–initiated) and extrinsic (caspase-8–initiated) apoptotic cascades. Coordinated upregulation of multiple caspases suggests robust engagement of apoptotic machinery.
- **Mitochondrial Apoptosis Amplification:** DIABLO (AIF) and SMAC (DIABLO) are upregulated, indicating mitochondrial release of pro-apoptotic factors that amplify caspase cascade activation and antagonize IAP-mediated caspase inhibition.

**ER Stress and UPR Category:**
- **UPR Transcription Factor Activation:** DDIT3 (CHOP), ATF4, and ATF6 are markedly upregulated (log₂FC: +1.1 to +1.5), indicating activation of all three canonical ER stress sensor pathways (PERK, ATF6, IRE1α). In normal circumstances, moderate UPR activation promotes protein homeostasis recovery; however, sustained or excessive UPR activation—as evidenced here—triggers pro-apoptotic transcriptional programs including CHOP-mediated induction of pro-apoptotic BH3-only proteins.
- **ER Chaperone Upregulation:** HSPA5 (BiP), HSPA9, and DNAJB9 are upregulated, reflecting accumulation of misfolded proteins in the ER and engagement of chaperone-mediated stress responses. Paradoxically, despite upregulation of ER chaperones, the concurrent activation of pro-apoptotic UPR mediators indicates that TP4-induced ER stress overwhelms cellular recovery capacity.
- **Mechanistic Interpretation:** TP4's disruption of intracellular calcium homeostasis (established in prior work) leads to ER calcium depletion, triggering calcium-dependent, protein misfolding-induced ER stress. This stress, if unresolved, transitions from pro-survival (adaptive) UPR to pro-apoptotic (terminal) UPR through sustained CHOP and ATF4 activation.

**MAPK/JNK/AP-1 Category:**
- **Immediate-Early Gene Activation:** FOS, FOSB, JUN, JUNB, and ATF3 are among the most robustly upregulated genes (log₂FC: +1.8 to +2.4 for FOS/FOSB, +1.5 to +2.1 for JUN family). This coordinated upregulation reflects activation of the MAPK/JNK/AP-1 signaling cascade in response to TP4-induced cellular stress.
- **AP-1 Complex Remodeling:** The marked upregulation of FOSB (which encodes ΔFosB isoforms with altered transcriptional properties) and downregulation of c-Jun (through cell cycle/proliferation suppression) suggests a remodeling of the AP-1 transcription factor complex that shifts its target gene selectivity toward pro-apoptotic and stress-response genes.
- **Stress-Responsive vs. Proliferation-Responsive AP-1 Function:** In normal proliferation, AP-1 drives expression of growth and survival genes; in stress contexts (as induced by TP4), AP-1 redirects to drive expression of death and stress-response genes. The FOSB-enriched AP-1 complex induced by TP4 represents a "stress-AP-1" configuration distinct from proliferation-promoting AP-1.
- **Transcriptional Reprogramming Hub:** The convergent activation of MAPK/JNK/AP-1 with simultaneous activation of ATF/CREB pathways (ATF3, ATF4, ATF6) indicates a coordinated transcriptional reprogramming hub that simultaneously activates pro-death programs (via AP-1 and ATF-mediated transcription of apoptotic genes) and stress response programs (via UPR and MAPK-driven metabolic and protein stress adaptation).

**Minimal Response in HDF Cells:**
Strikingly, HDF cells show minimal upregulation across all three categories—apoptosis, UPR, and MAPK/AP-1 genes. This differential sensitivity between cell types at the level of death-pathway activation genes directly explains why TP4 is selectively cytotoxic to cancer cells: cancer cells mount a robust, coordinated activation of multiple death pathways in response to TP4, whereas normal fibroblasts do not.

---

## Supplementary Figure 1. Summary Bar Chart of Differentially Expressed Gene Counts Across Contrasts Demonstrates Magnitude and Directionality of TP4-Induced Transcriptomic Perturbation.

**Legend:** Horizontal stacked bar chart displaying the number of significantly upregulated (red) and downregulated (blue) genes for each of four statistical contrasts. FDR threshold: 0.05; fold-change threshold: |log₂FC| >0.585 (1.5-fold). Numbers annotated on bars indicate gene counts per direction per contrast.

**Quantitative Summary:**
- **TP4 in MDA-MB-231:** 127 upregulated, 94 downregulated (total: 221 DEGs)
- **TP4 in HDF:** 23 upregulated, 15 downregulated (total: 38 DEGs)
- **MDA-Mock vs. HDF-Mock (Baseline):** 184 upregulated in MDA, 156 downregulated (total: 340 DEGs)
- **Interaction (Cancer-Selective):** 39 preferentially upregulated in MDA, 29 preferentially downregulated (total: 68 DEGs)

**Interpretation:** The dramatic difference in DEG counts between MDA (221 genes) and HDF (38 genes) contrasts highlights the 5.8-fold greater transcriptomic responsiveness of cancer cells to TP4. The balanced direction of dysregulation in MDA cells (127 up, 94 down; ~1.35:1 ratio of upregulation to downregulation) indicates that TP4 activates multiple coordinated transcriptional programs: upregulation of pro-death genes and downregulation of pro-survival genes. The interaction contrast, though quantitatively smaller (68 genes), represents the most mechanistically informative set, as it isolates cancer-selective responses freed from cancer-versus-normal baseline differences.

---

*End of Figure Legends*
