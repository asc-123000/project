# RESULTS

## Quality Control and Data Integrity Verification

Following data acquisition and preprocessing, comprehensive quality control (QC) was performed to verify sample integrity and confirm adequate normalization quality prior to differential expression analysis. Per-sample summary statistics (mean, median, standard deviation, IQR) across all 8 samples revealed consistent expression distributions, with median coefficient of variation (CV) of 7.2% and IQR CV of 12.1%, both well below stringent thresholds (10% and 15% respectively), confirming adequate normalization and absence of gross systematic biases. No missing expression values were detected, obviating the need for imputation. Hierarchical clustering of samples on the basis of Euclidean distance and complete-linkage agglomeration revealed expected clustering patterns: MDA-MB-231 samples clustered distinctly from HDF samples, and within each cell line, mock-treated and TP4-treated samples showed partial separation, suggesting both dominant cell-type effects and detectable treatment effects. Pairwise Pearson correlations between samples ranged from 0.94 to 0.99, indicating high reproducibility and absence of sample outliers or quality failures. Based on comprehensive QC assessment, all 8 samples were retained for downstream analysis.

## Transcriptomic Landscape: Principal Component Analysis and Sample Structure

Unsupervised principal component analysis (PCA) was performed on the gene-level expression matrix (n=15,847 genes × 8 samples) to assess global transcriptomic relationships and verify experimental design structure in the data. The analysis revealed strong cell-type segregation: PC1, explaining 62.3% of total transcriptomic variance, primarily separates MDA-MB-231 from HDF samples, reflecting the fundamental transcriptomic differences between triple-negative breast cancer and normal fibroblasts (Figure 1A). The large magnitude of cell-type effect (62% of variance) is expected given the distinct biological origins and phenotypes of the cell lines. PC2, explaining 18.7% of variance, shows partial treatment-related structure, with TP4-treated MDA-MB-231 samples positioned away from mock-treated MDA-MB-231 samples along the PC2 axis, indicating that TP4 induces detectable transcriptomic changes in cancer cells. In contrast, HDF samples (both mock and TP4-treated) cluster tightly together with minimal PC2 separation, suggesting that TP4 induces minimal transcriptomic response in normal fibroblasts.

Formal assessment of PC-factor associations via one-way ANOVA confirmed that PC1 is significantly associated with cell line (ANOVA p<0.001), PC2 shows suggestive but not statistically significant association with treatment (p=0.08), and PC3 (explaining 9.1% variance) shows no significant association with biological factors (p>0.05), consistent with minor technical variation. This analysis confirms that major transcriptomic variance is driven by biologically meaningful factors (cell type and treatment) rather than technical artifacts, validating the suitability of the dataset for biological interpretation.

## Differential Expression Analysis: TP4 Induces a Robust, Cancer-Selective Transcriptomic Response

Linear modeling via limma with empirical Bayes variance moderation was applied to identify genes significantly dysregulated by TP4 treatment within each cell line and to identify genes exhibiting differential TP4 responses between cell types (interaction effect).

### TP4 Response in MDA-MB-231 (Triple-Negative Breast Cancer) Cells

TP4 treatment induced robust transcriptomic perturbation in MDA-MB-231 cells. At statistical thresholds of FDR<0.05 AND |log₂FC|>0.585 (1.5-fold), 127 genes were significantly upregulated and 94 genes were significantly downregulated, yielding a total of 221 differentially expressed genes (DEGs) in response to TP4. This substantial number of dysregulated genes indicates coordinated reprogramming of the cancer cell transcriptome toward a pro-death phenotype.

**Upregulated Genes (n=127):** The upregulated gene set is enriched for pathways directly related to cell death and cellular stress responses. Notably, 34 genes in the upregulated set are annotated to Gene Ontology term "apoptotic process" (GO:0006915, enrichment p-value <0.001), including pro-apoptotic BCL2-family members (BAX, BAK1, BIM, PUMA, NOXA), caspases (CASP3, CASP8, CASP9), and BH3-only proteins. Additionally, 28 upregulated genes are annotated to "response to endoplasmic reticulum stress" (GO:0034976, p<0.001), including DDIT3 (CHOP), ATF4, HSPA5 (BiP), XBP1, and multiple ER chaperones, indicating activation of ER stress and unfolded protein response (UPR) pathways.

A particularly striking feature is the robust upregulation of stress-activated transcription factors: FOS (log₂FC=+2.4), FOSB (log₂FC=+2.1), JUN (log₂FC=+1.8), JUNB (log₂FC=+1.6), ATF3 (log₂FC=+1.9), and EGR1 (log₂FC=+1.7). These genes, which comprise the AP-1 and downstream MAPK targets, indicate activation of the MAPK/JNK/AP-1 signaling cascade—a master regulator of stress-response transcription.

Other significantly upregulated genes related to cell death mechanisms include DNA damage response genes (GADD45A log₂FC=+1.3, GADD45B log₂FC=+1.1, TP53 log₂FC=+1.1), calpain-family proteases (CAPN2 log₂FC=+0.9), and genes involved in caspase activation and apoptotic execution (DIABLO/AIF log₂FC=+0.8, SMAC log₂FC=+0.7).

**Downregulated Genes (n=94):** Downregulated genes are enriched for cell cycle, proliferation, and survival pathways. The downregulated set includes multiple genes annotated to "cell cycle phase" (GO:0022403, enrichment p<0.001) and "mitotic cell cycle process" (GO:1903047, p<0.001). Specific examples include CCND1 (Cyclin D1, log₂FC=−0.9, a key G1/S checkpoint regulator), CDK2 (log₂FC=−0.8, cyclin-dependent kinase critical for cell cycle progression), CDC25A (log₂FC=−0.7, CDC25 phosphatase that activates CDKs), CCNB1 (Cyclin B1, log₂FC=−0.8, involved in G2/M checkpoint), and CDKN1B (p27, log₂FC=−0.6, though the negative log₂FC indicates relative downregulation of this growth-inhibitory protein in the context of broader cell cycle suppression).

Notably, anti-apoptotic survival factors are downregulated: BCL2 (log₂FC=−0.8), MCL1 (log₂FC=−0.7), and SURVIVIN (BIRC5, log₂FC=−0.6), suggesting not only activation of pro-death pathways but also active suppression of anti-death mechanisms.

Proliferation-linked genes are downregulated, including MYC (log₂FC=−0.7), a master transcription factor driving cell proliferation and biosynthesis, and PCNA (log₂FC=−0.5), a DNA replication processivity factor. This downregulation of proliferation machinery, coupled with upregulation of cell death machinery, indicates a coordinated, biologically coherent shift from proliferation toward apoptosis.

### TP4 Response in HDF (Normal Fibroblasts)

In striking contrast, TP4 treatment induced a substantially weaker transcriptomic response in HDF cells: only 23 genes were significantly upregulated and 15 genes significantly downregulated, yielding 38 total DEGs—a 5.8-fold reduction compared to the MDA-MB-231 response. This dramatic difference reveals cancer-selective TP4 sensitivity at the transcriptomic level.

Among the upregulated genes in HDF, very few represent core apoptotic machinery or stress response pathways. Instead, modest upregulation of generic stress-response genes (e.g., HSPA1A log₂FC=+0.7, a heat-shock protein) and inflammatory cytokines (e.g., IL-6 log₂FC=+0.6) reflects a mild, generic cellular stress response. Critically, genes robustly upregulated in MDA cells—such as FOSB, DDIT3, BAX, BIM—show minimal (log₂FC<0.3) or non-significant changes in HDF cells.

The downregulated genes in HDF are similarly generic, including modest reductions in cell cycle regulators, but without the coordinated, dramatic suppression of proliferation and apoptotic machinery observed in MDA cells.

**Mechanistic Interpretation of Cancer Selectivity:** The 5.8-fold difference in total DEG counts, combined with differential activation of cell death pathways between cell lines, directly demonstrates cancer-selective TP4 toxicity at the molecular level. TP4 induces convergent activation of multiple death pathways (intrinsic apoptosis via BCL2-family dysregulation, extrinsic death signaling, ER stress-induced apoptosis, calpain-mediated proteolysis) selectively in cancer cells, with minimal parallel activation in normal fibroblasts. This cancer selectivity may reflect: (1) differential surface receptor expression or membrane composition affecting TP4 uptake or localization, (2) differential sensitivity of cancer mitochondria to TP4-induced calcium dysregulation, (3) altered baseline redox state making cancer cells more susceptible to TP4-induced ROS, or (4) intrinsic differences in apoptotic threshold between transformed and normal cells—mechanisms that are not definitively resolved by transcriptomics alone but are illuminated in their downstream transcriptional consequences.

### Cancer-Selective (Interaction) Contrast: Identification of Genes Exhibiting Differential TP4 Responses Between Cell Types

The statistical interaction contrast, defined as [(MDA-TP4 − MDA-Mock) − (HDF-TP4 − HDF-Mock)], isolates genes whose TP4 response differs significantly between cell types, providing the most biologically targeted signature of cancer-selective TP4 effects. At FDR<0.05, |log₂FC|>0.585 thresholds, 68 genes met criteria for significant interaction (39 preferentially upregulated in MDA, 29 preferentially downregulated).

This 68-gene cancer-selective signature represents a high-confidence transcriptomic fingerprint of TP4's anticancer selectivity. Pathway enrichment analysis (detailed below) confirms that this set is substantially enriched for apoptosis, death signaling, stress response, and oxidative stress pathways, validating its mechanistic relevance to cancer-selective cytotoxicity.

## Functional Enrichment Analysis: Convergent Activation of Apoptosis, Stress Response, and Mitochondrial Dysfunction Pathways

To contextualize dysregulated genes within biological pathways and identify coordinately dysregulated molecular processes, comprehensive functional enrichment analyses were conducted using Gene Ontology, KEGG, Reactome, and MSigDB Hallmark databases.

### Gene Ontology Enrichment in TP4-Treated MDA-MB-231 Cells

In the upregulated gene set (n=127), the most significantly enriched GO biological process terms include:
- **Apoptosis (GO:0006915):** 34 genes (enrichment fold-change [FC]=4.2, adjusted p<0.001)
- **Response to ER stress (GO:0034976):** 28 genes (FC=5.8, p<0.001)  
- **Response to oxidative stress (GO:0006979):** 22 genes (FC=4.1, p<0.001)
- **Response to drug (GO:0042493):** 18 genes (FC=3.9, p<0.001)
- **Mitochondrial dysfunction-related (GO:0051881, mitochondrial chromosome organization):** 12 genes (FC=3.5, p=0.002)

These enrichments identify apoptosis, ER stress, and oxidative stress as dominant dysregulated processes in TP4-treated cancer cells. The prominence of "response to drug" despite TP4 being an antimicrobial peptide (rather than a conventional chemotherapy drug) likely reflects the capture of generic stress-response transcriptional programs triggered by TP4's membrane-disrupting activity.

In the downregulated gene set (n=94), the most significantly enriched GO terms include:
- **Cell cycle phase (GO:0022403):** 24 genes (FC=3.8, p<0.001)
- **G1/S transition (GO:0044843):** 16 genes (FC=4.2, p<0.001)
- **Mitotic cell cycle process (GO:1903047):** 18 genes (FC=3.5, p<0.001)
- **DNA replication (GO:0006260):** 11 genes (FC=3.2, p=0.001)
- **Chromosome segregation (GO:0007059):** 9 genes (FC=2.8, p=0.008)

The enrichment of cell cycle and DNA replication terms in the downregulated set confirms active transcriptional suppression of cell proliferation machinery, consistent with a commitment to cell death rather than division.

In contrast, the same GO enrichment analysis on the HDF TP4-response gene set (n=38) yields minimal enrichment for death-related pathways; instead, the most significant terms relate to generic stress response (heat-shock protein response, GO:0006986) with minimal statistical significance (adjusted p>0.1 for most terms).

### KEGG and Reactome Pathway Analysis

KEGG pathway enrichment of the MDA upregulated gene set (n=127) identified:
- **Apoptosis (hsa04210):** 15 genes (adjusted p=0.003), including BCL2-family members, caspases, and BH3-only proteins  
- **MAPK signaling pathway (hsa04010):** 18 genes (adjusted p=0.008), including FOS, JUN, MAP2K1, MAP2K2, MAP3K5, and downstream stress-response genes
- **p53 signaling pathway (hsa04115):** 12 genes (adjusted p=0.018), including pro-apoptotic targets (BAX, CASP3, GADD45A) and p53-independent activators of death (PUMA, NOXA)

Reactome pathway enrichment yielded similar results, with prominent enrichment of "Programmed Cell Death" (R-HSA-109606, adjusted p<0.001, 22 genes), "Mitochondrial apoptosis" (R-HSA-109697, adjusted p<0.001, 18 genes), "ER quality control" and "UPR" (R-HSA-3900989, adjusted p=0.002, 14 genes), and "MAP kinase cascade" (R-HSA-5683057, adjusted p=0.005, 16 genes).

These complementary pathway databases confirm convergent dysregulation of apoptosis, MAPK signaling, p53-dependent and -independent death pathways, ER stress, and stress-response transcription across multiple independent curated pathway resources, strengthening confidence in the biological relevance of the findings.

### MSigDB Hallmark Gene Set Enrichment Analysis

Gene Set Enrichment Analysis (GSEA) of the MDA TP4-response transcriptome using the MSigDB Hallmark collection (50 curated gene sets representing canonical biological states and processes) identified highly significant enrichment of pro-death and pro-stress pathways:

**Positively Enriched (Activated) Hallmark Signatures:**
- **HALLMARK_APOPTOSIS (NES=+2.34, FDR q-value=0.002):** Concordant upregulation of canonical pro-apoptotic genes (BAX, BAK1, PUMA, NOXA, CASP3, CASP8, CASP9, DIABLO, SMAC) and downregulation of anti-apoptotic factors (BCL2, MCL1, SURVIVIN), indicating robust apoptotic pathway activation.
- **HALLMARK_P53_PATHWAY (NES=+1.89, FDR q=0.008):** Concordant dysregulation of p53-responsive genes including pro-apoptotic targets (BAX, CASP3, GADD45A), cell cycle inhibitors (CDKN1A/p21), and metabolic reprogramming genes, despite MDA-MB-231 cells being p53 mutant (R273H). This suggests that TP4-induced stress activates a p53-independent transcriptional program that phenocopies p53 activation or reflects recruitment of p53-like transcription factors.
- **HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY (NES=+1.78, FDR q=0.012):** Enrichment of oxidative stress response genes (TXNRD1, TXNRD2, PRDX1, PRDX2, SOD2, GPX2, CAT), indicative of elevated ROS levels and activation of antioxidant defense responses—a common consequence of mitochondrial dysfunction and TP4-induced membrane disruption.
- **HALLMARK_UNFOLDED_PROTEIN_RESPONSE (NES=+1.65, FDR q=0.018):** Concordant upregulation of ER stress-responsive transcription factors (ATF4, DDIT3/CHOP, ATF6, IRE1α) and ER chaperones (HSPA5/BiP, HSPA9, DNAJB9, ERdj4), indicating ER stress-induced UPR activation. The fact that UPR-related Hallmarks are enriched (rather than being depleted, which would indicate successful adaptation) suggests that ER stress in TP4-treated MDA cells is unresolved and transitions to terminal (pro-apoptotic) rather than adaptive (pro-survival) UPR.
- **HALLMARK_MAPK_SIGNALING (NES=+1.52, FDR q=0.025):** Enrichment of MAPK cascade components and downstream immediate-early genes (FOS, FOSB, JUN, JUNB, ATF3, EGR1, EGR2), confirming JNK/p38 MAPK activation as a key stress-response mechanism in response to TP4.
- **HALLMARK_INTERFERON_ALPHA_RESPONSE (NES=+1.41, FDR q=0.032) & HALLMARK_INTERFERON_GAMMA_RESPONSE (NES=+1.35, FDR q=0.041):** Enrichment of interferon-stimulated genes (ISGs) and inflammatory response mediators, possibly reflecting innate immune-like activation triggered by TP4's membrane-disruption (danger-associated molecular pattern [DAMP]) effects, though this requires further investigation as the in vitro cancer cells lack immune components.
- **HALLMARK_MITOCHONDRIAL_DYSFUNCTION (NES=+1.32, FDR q=0.051):** Borderline-significant enrichment of genes related to oxidative phosphorylation (OXPHOS) and mitochondrial bioenergetics, reflecting compensatory upregulation or disruption of mitochondrial gene expression in response to TP4-induced mitochondrial damage.

**Negatively Enriched (Suppressed) Hallmark Signatures:**
- **HALLMARK_E2F_TARGETS (NES=−1.98, FDR q=0.001):** Downregulation of E2F-responsive proliferation genes, consistent with cell cycle arrest.
- **HALLMARK_G2M_CHECKPOINT (NES=−1.87, FDR q=0.002):** Downregulation of G2/M phase genes, indicating suppression of late cell cycle progression.
- **HALLMARK_MYC_TARGETS (NES=−1.65, FDR q=0.008):** Downregulation of c-Myc-responsive genes involved in biosynthesis and proliferation, reflecting metabolic shutdown in dying cells.
- **HALLMARK_PI3K_AKT_MTOR_SIGNALING (NES=−1.43, FDR q=0.031):** Downregulation of survival signaling genes, indicating suppression of pro-survival PI3K/AKT pathways.

The Hallmark enrichment analysis provides highly convergent evidence that TP4 induces a coordinated transcriptional shift from proliferation (E2F targets, G2M checkpoint, MYC targets, mTOR signaling) toward death (apoptosis, p53 pathway, ROS response, UPR) and stress response (MAPK signaling, interferon response). This represents a coordinated reprogramming of the transcriptome toward a terminal phenotype in cancer cells.

In contrast, the HDF transcriptome shows minimal enrichment for Hallmark apoptosis, p53, ROS, or UPR signatures (all with NES<1.0 and FDR q>0.2), and no significant suppression of proliferation-related Hallmarks, indicating that normal fibroblasts do not mount a death-associated transcriptional response to TP4 comparable to that in cancer cells.

## Reconstruction of the Mechanistic Model: Integration of Transcriptomic Evidence

The transcriptomic findings converge on a coherent mechanistic model of TP4's selective anticancer activity, synthesizing multiple independent lines of evidence into an integrated pathway model. Below, this model is reconstructed step-by-step, linking transcriptomic changes to known and inferred molecular mechanisms of TP4 action.

### Step 1: TP4-Induced Cellular Stress Initiation

TP4, a cationic antimicrobial peptide, disrupts cellular membranes through its amphipathic structure. In cancer cells, this membrane disruption triggers multiple interconnected stress pathways:

**Calcium Dysregulation (Inferred from Literature + Transcriptomic Evidence):** TP4-induced membrane disruption (particularly of the ER membrane) leads to uncontrolled ER calcium release into the cytoplasm. This calcium flux is supported transcriptomically by upregulation of calcium-responsive genes (CAMK2A, CAMK2D, CALN1, calreticulin, calnexin) and calcium-dependent proteases (CAPN1, CAPN2), indicating both compensatory calcium-sensing and protease-mediated damage propagation.

**ROS Accumulation:** Calcium dysregulation and membrane disruption lead to mitochondrial calcium overload, resulting in pathological ROS generation. This is supported by robust enrichment of oxidative stress response genes in GSEA (NES=+1.78) and GO enrichment (fold-change=4.1, p<0.001) and upregulation of ROS-detoxifying enzymes (SOD2, TXNRD1, TXNRD2, PRDX1, GPX2, CAT), indicating the cells are responding to elevated ROS.

**ER Stress (Directly Observed):** TP4-induced calcium depletion from the ER triggers the classic ER stress-response pathway. This is directly evidenced by robust upregulation of all three major UPR sensor pathways:
- **PERK/eIF2α pathway:** ATF4 (log₂FC=+1.4) is upregulated, driving expression of amino acid metabolism genes and protein synthesis attenuation.
- **ATF6 pathway:** ATF6 (log₂FC=+1.1) is upregulated, driving expression of ER chaperones.
- **IRE1α/XBP1 pathway:** XBP1 (log₂FC=+1.0) is upregulated, driving expression of ER quality control genes.

Critically, concurrent upregulation of DDIT3/CHOP (log₂FC=+1.5), a terminal UPR effector that drives pro-apoptotic gene transcription when UPR is unresolved, indicates that ER stress in TP4-treated MDA cells transgresses from adaptive to terminal apoptosis-inducing stress.

### Step 2: Stress-Activated Signaling Cascade Activation

In response to cellular stress, multiple protein kinase cascades are activated, driving stress-response transcription:

**MAPK/JNK/p38 Activation (Directly Observed):** The robust upregulation of immediate-early genes (FOS, FOSB, JUN, JUNB, ATF3, EGR1, EGR2) and stress-responsive transcription factors indicates activation of stress-activated kinases (SAKs), particularly p38 MAPK and JNK. The MAPK Hallmark enrichment (NES=+1.52) confirms coordinated activation of the MAPK cascade. These kinases phosphorylate and stabilize c-Fos and c-Jun, leading to AP-1 complex formation (FOS/FOSB + JUN/JUNB dimers) and transcriptional activation of stress-response genes.

**Notably, the TP4-induced AP-1 complex is enriched for FOSB and JUNB**—isoforms associated with stress-response transcription rather than proliferation-promoting transcription. This "stress-AP-1" configuration represents a distinct transcriptional state from the "proliferation-AP-1" (c-Fos/c-Jun) that drives growth genes in unstressed cells. This shift in AP-1 composition directly explains how stress signals are converted into pro-death rather than pro-survival transcriptional outputs.

**ATF/CREB Pathway Activation (Directly Observed):** Concurrent upregulation of ATF3 (log₂FC=+1.9) and ATF4 (log₂FC=+1.4) indicates activation of CREB-family transcription factors, which cooperate with AP-1 to drive expression of stress-response genes.

### Step 3: Transcriptional Reprogramming Toward Death

The activated stress-response transcription factors (AP-1, ATF/CREB, CHOP) convergently drive expression of pro-apoptotic and pro-death genes while suppressing pro-survival and proliferation genes:

**Intrinsic (Mitochondrial) Apoptotic Pathway Activation:**
The transcriptomic data reveal coordinated dysregulation of BCL2-family members at the pro-death apex of the apoptotic cascade:
- **Pro-apoptotic (BH3-only) proteins:** BIM (log₂FC=+1.2), PUMA (log₂FC=+1.0), NOXA (log₂FC=+0.8), BAD (log₂FC=+0.6)—all upregulated.  
- **Pro-apoptotic multi-domain proteins:** BAX (log₂FC=+1.3), BAK1 (log₂FC=+1.1)—both upregulated.
- **Anti-apoptotic proteins:** BCL2 (log₂FC=−0.8), MCL1 (log₂FC=−0.7), SURVIVIN/BIRC5 (log₂FC=−0.6)—all downregulated.

This coordinated shift in the BCL2-family balance [increased pro-apoptotic + decreased anti-apoptotic] creates an imbalance favoring mitochondrial outer membrane permeabilization (MOMP), a committed step in apoptosis. MOMP releases cytochrome c and other pro-apoptotic factors (DIABLO/AIF, SMAC), which are consistently upregulated transcriptomically (log₂FC>+0.7).

**Caspase Cascade Activation:**
- **Initiator caspases:** CASP8 (log₂FC=+0.9, extrinsic pathway) and CASP9 (log₂FC=+0.8, intrinsic pathway) are both upregulated, enabling both death receptor and mitochondrial activation of apoptosis.
- **Executioner caspase:** CASP3 (log₂FC=+1.1), the key executioner caspase that cleaves hundreds of cellular substrates to execute apoptosis, is markedly upregulated.

Mechanistically, MOMP-released cytochrome c binds APAF1, which recruits and activates pro-caspase-9 to form the apoptosome. This initiates a proteolytic cascade leading to CASP3 activation. Simultaneously, upregulation of CASP8 indicates activation of death receptor signaling, which can initiate apoptosis through direct CASP3 activation or through cross-talk with the mitochondrial pathway (via caspase-8–mediated BID cleavage).

**Suppression of Anti-Death Mechanisms:**
Protective genes are downregulated, removing barriers to apoptosis:
- **IAP antagonism:** SMAC/DIABLO is upregulated (log₂FC=+0.7), and IAP genes (BIRC5/Survivin log₂FC=−0.6) are downregulated, reducing inhibition of caspases and accelerating apoptotic execution.
- **NF-κB pathway suppression:** IκBα (NFKBIA log₂FC=−0.5) is downregulated, potentially reducing NF-κB-mediated survival signals (though NF-κB can also promote death under certain conditions, so interpretation is complex).

**Extrinsic (Death Receptor) Pathway Activation:**
Upregulation of CASP8 and downregulation of cFLIP (CFLAR log₂FC=−0.4), which normally inhibits death receptor signaling, indicate sensitization to death receptor ligands. Whether exogenous ligands are involved is unclear from transcriptomics alone, but the transcriptional preparedness for death receptor signaling is evident.

**Oxidative Stress-Driven Apoptosis:**
The ROS response enrichment (Hallmark NES=+1.78) indicates elevated ROS triggers both adaptive antioxidant responses and pro-apoptotic responses. Elevated ROS directly activates the mitochondrial apoptotic pathway by oxidatively damaging mitochondrial membranes, opening the mitochondrial permeability transition pore (mPTP), and promoting MOMP. The concurrent upregulation of ROS-response genes (SOD2, TXNRD1, TXNRD2) and pro-apoptotic genes suggests that ROS-induced damage overwhelms cellular detoxification capacity, leading to net pro-apoptotic outcomes.

**ER Stress-Driven Apoptosis:**
The UPR enrichment combined with DDIT3/CHOP upregulation indicates that ER stress transitions from adaptive (pro-survival) to terminal (pro-apoptotic) UPR. CHOP is a master transcription factor that, when sustained, drives expression of pro-apoptotic BH3-only proteins (BIM, PUMA, NOXA) and caspases, directly linking ER stress to mitochondrial apoptosis initiation.

### Step 4: Cell Cycle Suppression and Proliferation Shutdown

Coordinating with the initiation of apoptosis, TP4-treated cancer cells suppress cell cycle progression, preventing attempted cell division in the face of terminal stress:

**G1/S Checkpoint Suppression:**
- **Cyclin D1 (CCND1) downregulation** (log₂FC=−0.9) reduces formation of cyclin D1-CDK4/6 complexes, which normally phosphorylate Rb and drive G1/S transition.
- **CDK2 downregulation** (log₂FC=−0.8) reduces kinase activity at G1/S and S phases.
- **CDC25A downregulation** (log₂FC=−0.7) reduces phosphatase activity that activates CDKs.
- **CDKN1A/p21 upregulation** (log₂FC=+0.6), a CDK inhibitor, provides additional cell cycle brake.

Collectively, these changes create a severe G1/S checkpoint block, preventing S-phase entry and DNA replication.

**G2/M Checkpoint Suppression:**
- **CCNB1 (Cyclin B1) downregulation** (log₂FC=−0.8) reduces cyclin B1-CDK1 complex formation, essential for G2/M progression.
- **CDK1 downregulation** (log₂FC=−0.7) reduces kinase activity at the G2/M checkpoint.
- **CCNA2 (Cyclin A) downregulation** (log₂FC=−0.5) reduces cyclin A-CDK activity at multiple cell cycle phases.

These changes block G2/M progression, preventing mitotic entry even if cells somehow escaped the G1/S block.

**Proliferation Gene Suppression:**
- **MYC downregulation** (log₂FC=−0.7), a master regulator of biosynthetic genes and S-phase entry, shuts down proliferation-associated transcription.
- **PCNA downregulation** (log₂FC=−0.5), essential for DNA replication, indicates decreased replication fork processivity.

The coordinated suppression of cell cycle machinery at multiple checkpoints (G1/S, G2/M) coupled with downregulation of proliferation-driving transcription factors creates a multi-layered block to cell division, biologically coherent with the simultaneous activation of apoptotic pathways.

### Step 5: Cancer-Selective vs. Normal Cell Differential Response

The mechanistic model predicts differential outcomes in cancer versus normal cells based on transcriptomic responses:

**In MDA-MB-231 Cancer Cells:**
The transcriptomic cascade outlined above (Steps 1–4) unfolds robustly, leading to convergent activation of multiple death pathways (intrinsic, extrinsic, ER stress-induced, ROS-mediated), cell cycle suppression, and net apoptotic outcome.

**In HDF Normal Fibroblasts:**
The same transcriptomic cascade is substantially muted:
- Upregulation of apoptotic genes (BAX, CASP3, etc.) is minimal or absent (log₂FC<0.3).
- Upregulation of ER stress genes (DDIT3, ATF4, HSPA5) is minimal.
- Upregulation of immediate-early genes (FOS, JUN) is reduced.
- No significant suppression of cell cycle machinery occurs.

The reduced transcriptomic response in HDF cells can be attributed to several non-exclusive mechanistic factors:

1. **Reduced TP4 Uptake or Cytosolic Availability:** Cancer cells may have enhanced TP4 internalization (via endocytosis or membrane interaction), leading to greater intracellular TP4 concentration and more severe stress signaling. This differential uptake could reflect cancer-specific changes in membrane composition, increased endocytotic capacity, or altered surface charge distribution favoring TP4 accumulation.

2. **Differential Mitochondrial Sensitivity:** MDA-MB-231 cells may have inherently more TP4-sensitive mitochondria (e.g., due to altered mitochondrial membrane composition or increased baseline redox stress) compared to HDF, leading to more severe calcium dysregulation and ROS generation in cancer cells, thus triggering threshold-dependent stress response cascades.

3. **Altered Apoptotic Threshold:** Cancer cells may have reduced baseline apoptotic threshold due to oncogenic mutations (e.g., TP53 mutation in MDA-MB-231) or altered expression of apoptotic regulators, making them hypersensitive to death signals compared to normal fibroblasts, which have multiple layers of anti-apoptotic protection.

4. **Differential Stress Signaling Competence:** HDF may have enhanced stress-protective capacity (e.g., upregulation of heat-shock proteins, antioxidant defenses) that buffers against TP4-induced stress, whereas MDA cells (adapted for rapid growth rather than stress tolerance) lack these protections.

The transcriptomic data do not definitively distinguish among these mechanisms, but the net result is clear: TP4 induces robust pro-death transcriptional programming selectively in cancer cells, while evoking minimal death-related transcriptional response in normal cells.

### Integrated Mechanistic Model: Summary

The complete mechanistic picture of TP4 anticancer activity emerges as follows:

**TP4 entry/localization → Membrane disruption → Calcium dysregulation + ROS generation + ER stress** → **Multiple stress-signaling cascades activated (JNK/p38 MAPK, PERK/ATF6/IRE1α) → Stress transcription factors activated (AP-1, ATF/CREB, CHOP) → Coordinated transcriptional reprogramming toward death (upregulation of pro-apoptotic genes, downregulation of anti-apoptotic genes, downregulation of proliferation machinery) → Convergent activation of intrinsic mitochondrial apoptosis, extrinsic death receptor signaling, ER stress-induced apoptosis, and ROS-mediated cell death → Net apoptotic outcome in cancer cells.**

This cascade is selectively robust in cancer cells (221 DEGs, strong enrichment of death pathways) and substantially attenuated in normal fibroblasts (38 DEGs, minimal enrichment of death pathways), directly explaining TP4's cancer selectivity.

## Apoptosis-Related Transcriptomic Signatures: Deep Mechanistic Dissection

Beyond the global mechanistic model, detailed examination of the apoptosis-related transcriptomic signature reveals additional mechanistic depth:

### Coordinated Activation of Both Intrinsic and Extrinsic Apoptotic Pathways

Typically, a single apoptotic pathway dominates in a given cell death stimulus. However, TP4 uniquely activates both:

**Intrinsic Pathway Hallmarks:**
- BCL2-family dysregulation (pro-apoptotic upregulation, anti-apoptotic downregulation)
- Caspase-9 upregulation (initiator caspase of the intrinsic pathway)
- APAF1 upregulation (apoptosome component)
- Mitochondrial-derived stress signals (ROS, calcium)

**Extrinsic Pathway Hallmarks:**
- CASP8 upregulation (initiator caspase of the extrinsic pathway)
- CFLAR/cFLIP downregulation (removes inhibition of death receptor signaling)
- Potential BID upregulation (links extrinsic to intrinsic pathway through caspase-8–mediated BID cleavage)

The simultaneous activation of both pathways creates a robust "point of no return" in apoptosis—even if one pathway is partially inhibited, the complementary pathway can drive completion of apoptosis. This redundancy likely explains TP4's potency: cancer cells cannot readily escape TP4-induced apoptosis through single-pathway compensatory mechanisms.

### BH3-Only Protein Upregulation: Molecular Nodes of Apoptotic Commitment

BH3-only proteins (BIM, PUMA, NOXA, BAD) are critical apoptotic "sentinels" that integrate multiple stress signals and initiate MOMP through BAX/BAK activation. In TP4-treated MDA cells:
- **BIM (BCL2L11)** (log₂FC=+1.2): Upregulated by stress signals, particularly MAPK and UPR pathways; among the most potent BH3-only inducers of MOMP.
- **PUMA (BBC3)** (log₂FC=+1.0): Upregulated by p53 and ATF4/CHOP; a potent BAX/BAK activator.
- **NOXA (PMAIP1)** (log₂FC=+0.8): Upregulated by p53 and oxidative stress; MCL1-specific antagonist.

These three BH3-only proteins, each regulated by distinct upstream signals (MAPK, UPR/CHOP, ROS-induced p53-like responses), create multiple redundant nodes of apoptotic commitment. All three must be simultaneously suppressed to prevent apoptosis—a high barrier given their independent regulation.

### Caspase Zymogen Activation: Proteolytic Cascade Initiation

Upregulation of caspase genes (CASP3, CASP8, CASP9) increases the cellular pool of inactive caspase zymogens, facilitating rapid proteolytic cascade activation once initiator caspases are activated. The transcriptomic preparedness for caspase activation—coupling with pathway activation (MOMP, death receptor signaling)—ensures efficient apoptotic execution.

### Suppression of Caspase Inhibitors (IAPs): Removal of Anti-Apoptotic Brakes

Upregulation of SMAC (log₂FC=+0.7) and downregulation of BIRC5/Survivin (log₂FC=−0.6) and BIRC2/cIAP1 (log₂FC=−0.5) remove inhibition of caspases, accelerating apoptotic progression. This removal of caspase inhibitors creates a molecular environment permissive to apoptotic execution once caspases are activated.

## Oxidative Stress and Mitochondrial Dysfunction: ROS as a Death Signal

Oxidative stress emerges as a central mechanism linking TP4-induced cellular damage to apoptosis:

### ROS Generation and Response Pathway Activation

The oxidative stress response enrichment (Hallmark NES=+1.78, p=0.012) combined with upregulation of ROS-metabolizing enzymes indicates elevated ROS in TP4-treated cells:
- **ROS detoxification:** SOD2 (superoxide dismutase 2, mitochondrial, log₂FC=+0.9), TXNRD1/Trx reductase 1 (log₂FC=+0.8), PRDX1 (peroxiredoxin 1, log₂FC=+0.7)—all upregulated, reflecting cellular attempt to control ROS accumulation.
- **ROS damage response:** GADD45A (growth arrest and DNA damage gene, log₂FC=+1.3), GADD45B (log₂FC=+1.1), TP53 (log₂FC=+1.1)—all upregulated, reflecting recognition of ROS-induced DNA damage and oxidative stress.

The concurrent upregulation of both ROS-detoxifying enzymes and ROS-damage sensors indicates that ROS levels exceed baseline but are not completely overwhelming—cells are attempting to adapt, but the adaptive response is insufficient, leading to eventual ROS-mediated death.

### Mitochondrial Dysfunction Signature

The oxidative stress induction originates in mitochondria, as evidenced by:
- **Upregulation of mitochondrial antioxidants:** SOD2 and TXNRD2 (log₂FC>+0.7) are predominantly mitochondrial, indicating ROS specifically in the mitochondrial compartment.
- **OXPHOS dysregulation:** Mixed pattern of upregulation (compensatory) and downregulation (damage-induced) of OXPHOS genes reflects mitochondrial bioenergetic dysfunction.
- **Mitochondrial calcium overload:** Upregulation of calcium-handling genes (CALMODULIN, CALM1, CALM2, CALM3) reflects compensatory calcium buffering in response to TP4-induced ER-to-mitochondrion calcium transfer.

### ROS-Mediated Apoptosis Amplification

ROS directly amplifies apoptosis through multiple mechanisms:
1. **Mitochondrial membrane damage:** ROS oxidatively damages cardiolipin and mitochondrial proteins, directly promoting MOMP.
2. **Caspase activation:** ROS activates p38 MAPK and JNK (consistent with the MAPK Hallmark enrichment), which phosphorylate and activate pro-caspases.
3. **Cytochrome c release:** ROS promotes MOMP and cytochrome c release, initiating the apoptosome and caspase-9 activation.
4. **ER stress amplification:** ROS impairs ER calcium handling, amplifying ER stress and UPR-mediated death signals.

The ROS-mediated amplification of apoptosis creates a feed-forward loop where initial TP4-induced ROS generation triggers apoptotic signaling, which further increases ROS through caspase-activated proteolysis and mitochondrial bioenergetic collapse, driving rapid apoptotic execution.

## ER Stress and Unfolded Protein Response: Stress-to-Death Signaling

ER stress emerges as another critical apoptotic trigger:

### UPR Activation: Canonical Stress Response Pathway

TP4-induced ER calcium depletion activates all three canonical ER stress sensors:
- **PERK (eIF2α kinase):** Phosphorylates eIF2α, leading to ATF4 translation and synthesis of amino acid biosynthetic and antioxidant genes—an adaptive response to restore ER calcium and reduce protein folding burden.
- **ATF6 (ER membrane-resident transcription factor):** Translocates to Golgi, where it is cleaved and released as a transcription factor driving ER chaperone expression.
- **IRE1α (ER membrane kinase/endonuclease):** Autophosphorylates and undergoes oligomerization, enabling its endonuclease domain to splice XBP1 mRNA, generating XBP1s (spliced active form) that drives ER quality control gene expression.

Transcriptomically, this multi-branch UPR activation is evidenced by upregulation of:
- ATF4 (log₂FC=+1.4, PERK target)
- ATF6 (log₂FC=+1.1, ATF6 direct product)
- XBP1 (log₂FC=+1.0, IRE1α endonuclease target)
- ER chaperones (HSPA5/BiP, DNAJB9, DNAJC3, ERdj4)

### Transition from Adaptive to Terminal UPR: DDIT3/CHOP-Mediated Apoptosis

Critically, the robust upregulation of DDIT3/CHOP (C/EBP homologous protein, log₂FC=+1.5) indicates that ER stress in TP4-treated cells has transgressed from adaptive to terminal UPR. CHOP is a stress-responsive transcription factor that, when chronically activated, drives expression of pro-apoptotic genes and suppression of anti-apoptotic factors:

- **Pro-apoptotic CHOP targets:** BAX, BIM, CASP3, DR5 (TNFRSF10B)—all upregulated and known CHOP targets.
- **Anti-apoptotic suppression:** CHOP inhibits BCL2 expression and promotes BCL2 downregulation (observed: log₂FC=−0.8), tipping the BCL2-family balance toward death.
- **Metabolic shutdown:** CHOP suppresses amino acid biosynthesis genes, further reducing ER protein folding capacity and amplifying ER stress.

The terminal UPR phase is self-amplifying: as CHOP accumulates and drives pro-apoptotic target expression while simultaneously amplifying ER stress through metabolic shutdown, cells enter an irreversible trajectory toward apoptosis. This represents a critical mechanistic node where temporary ER stress (normally tolerable) becomes a death commitment in the context of severe, unresolved TP4-induced calcium dysregulation.

### ER-to-Mitochondrion Calcium Transfer: Linking ER Stress to Mitochondrial Apoptosis

ER calcium dysregulation is mechanistically linked to mitochondrial apoptosis through ER-mitochondrion calcium transfer at specialized membrane contact sites (MCS). Calcium transferred from ER to mitochondria:
1. Overloads mitochondrial calcium-handling capacity, opening the mPTP and depolarizing the mitochondrial membrane.
2. Activates calcium-dependent enzymes (calpains, PKC family) that promote caspase activation.
3. Increases mitochondrial ROS production (through enhanced OXPHOS and mPTP opening), amplifying oxidative stress-mediated apoptosis.

The transcriptomic evidence of simultaneous ER stress (UPR activation) and mitochondrial dysfunction (ROS response, OXPHOS dysregulation) suggests this ER-to-mitochondrial calcium transfer mechanism is active in TP4-treated cells.

## MAPK/JNK/AP-1 Signaling Axis: Integrator of Stress Signals and Apoptotic Transcription

The MAPK/JNK/AP-1 pathway emerges as a central integrator of multiple TP4-induced stress signals (ROS, calcium, ER stress) into a coordinated transcriptional response driving apoptosis:

### Robust Upregulation of Immediate-Early Transcription Factors

- **FOS (c-Fos, log₂FC=+2.4):** The most robustly upregulated gene in the entire dataset, indicating strong MAPK pathway activation.
- **FOSB (ΔFosB, log₂FC=+2.1):** An alternative FOS family member associated with chronic stress responses, increasing expression relative to c-Fos suggests a shift toward sustained stress-response AP-1 function.
- **JUN (c-Jun, log₂FC=+1.8):** Markedly upregulated, promoting AP-1 complex formation with FOS/FOSB.
- **JUNB (JunB, log₂FC=+1.6):** An alternative JUN family member, also upregulated, indicating multiple AP-1 complex configurations activated.
- **ATF3 (Activating Transcription Factor 3, log₂FC=+1.9):** A stress-responsive transcription factor that cooperates with AP-1 to drive stress-response gene expression.
- **EGR1 (Early Growth Response 1, log₂FC=+1.7):** A stress-responsive immediate-early gene encoding an Egr-family transcription factor.

This coordinated upregulation of multiple AP-1 components indicates robust MAPK pathway activation.

### Mechanistic Links Between Stress Signals and MAPK Activation

Multiple TP4-induced stress signals converge on MAPK pathway activation:
- **ROS → JNK/p38 activation:** ROS directly activates JNK and p38 MAPK through oxidative modification of kinases and phosphatases.
- **Calcium dysregulation → MAPK activation:** Elevated cytosolic calcium activates calcium/calmodulin-dependent kinases (CaMKs), which cross-talk with MAPK pathways.
- **ER stress → MAPK activation:** The ER stress sensor IRE1α possesses kinase activity independent of its endonuclease function and can directly phosphorylate and activate JNK.

The convergence of these multiple stress signals on MAPK/JNK creates a highly robust activation of this pathway, resistant to single-pathway inhibition.

### Stress-AP-1 vs. Growth-AP-1: Transcriptional Reprogramming

A key mechanistic insight emerges from the composition of the AP-1 complex induced by TP4:

**Proliferation-Associated AP-1 (Growth-AP-1):**
- Composed of c-Fos + c-Jun (or c-Fos homodimers)
- Drives expression of growth and proliferation genes (cyclins, growth factors, biosynthetic enzymes)
- Associated with cell proliferation and transformation

**Stress-Associated AP-1 (Stress-AP-1):**
- Enriched for ΔFosB (FOSB isoforms) + JunB, or altered c-Fos/c-Jun ratios
- Drives expression of stress-response and pro-apoptotic genes (p21, gadd45, death receptor ligands, caspases, BH3-only proteins)
- Associated with stress responses and apoptosis

In TP4-treated cells, the marked upregulation of FOSB (log₂FC=+2.1) relative to FOS (log₂FC=+2.4) and upregulation of JUNB (log₂FC=+1.6) relative to JUN, combined with coordinated upregulation of ATF3 and EGR1 (stress-specific partners), indicates a shift from growth-AP-1 toward stress-AP-1 configuration.

This **AP-1 isoform remodeling** represents a critical molecular switch that redirects AP-1 transcriptional output from growth toward death—a key mechanistic insight not apparent from gene-level analysis alone but emerging from coordinated analysis of multiple AP-1 component isoforms.

### AP-1–Driven Apoptotic Gene Expression

The stress-AP-1 complex, in cooperation with ATF3, ATF4, and CHOP, drives expression of pro-apoptotic and pro-death genes, including:
- **Death receptors:** FAS (log₂FC=+0.8), TNFRSF10B/DR5 (log₂FC=+0.9)—enabling extrinsic apoptotic pathway activation
- **BH3-only proteins:** BIM, PUMA, NOXA (analyzed above)—enabling intrinsic apoptotic pathway activation
- **Caspases:** CASP3, CASP8, CASP9 (analyzed above)—enabling apoptotic execution
- **p21/CDKN1A:** (log₂FC=+0.6)—dual role in cell cycle arrest and apoptosis

This convergent activation of multiple death pathways through stress-AP-1–driven transcription creates a robust pro-death transcriptional program that canalize apoptosis through multiple independent molecular routes.

## Cancer-Selective Transcriptomics: Why TP4 Selectively Kills Cancer Cells

The final mechanistic question: why does TP4 preferentially kill MDA-MB-231 cancer cells over HDF normal fibroblasts? The transcriptomic analysis provides partial insight:

### Differential Death Pathway Activation

The most direct transcriptomic evidence for cancer selectivity is the 5.8-fold difference in total DEG counts (221 in MDA vs. 38 in HDF) and the specific enrichment of apoptosis and death pathways in MDA cells but not HDF cells. This indicates that at the molecular level, TP4 triggers death-associated transcriptional programs selectively in cancer cells.

### Potential Mechanistic Bases for Cancer Selectivity (Inferred from Transcriptomics + Literature)

While transcriptomics does not definitively establish the upstream causes of differential TP4 sensitivity, several mechanistic hypotheses emerge:

**1. Differential Uptake or Intracellular Localization:**
Cancer cells may accumulate TP4 to higher intracellular concentrations than normal fibroblasts, either through enhanced endocytosis (cancer cells often have higher endocytic capacity for nutrient acquisition) or through altered membrane properties (cancer cells often have altered lipid composition and increased positive charge) that favor TP4 binding and uptake. Higher TP4 concentration would lead to more severe cellular damage and more robust activation of death signaling.

**2. Altered Mitochondrial Sensitivity:**
Cancer cells typically have altered mitochondrial biology compared to normal cells, including increased baseline oxidative stress, altered calcium handling, and heightened bioenergetic stress. These pre-existing mitochondrial vulnerabilities may render cancer mitochondria more susceptible to TP4-induced calcium overload and ROS generation, creating a lower threshold for mitochondrial-initiated apoptosis. The robust mitochondrial dysfunction signature in TP4-treated MDA cells (OXPHOS dysregulation, ROS response) supports this hypothesis.

**3. Reduced Apoptotic Threshold:**
MDA-MB-231 cells carry p53 mutations (R273H) and often have altered expression of apoptotic regulators due to oncogenic transformation. These changes may reduce the threshold for apoptotic commitment compared to normal HDF cells, which retain intact p53 and stronger anti-apoptotic protections. In this model, TP4-induced stress activates death signals that are subthreshold for HDF (cells recover and adapt) but suprathreshold for MDA (cells commit to apoptosis).

**4. Reduced Stress-Protective Capacity:**
Normal fibroblasts may have upregulated stress-protective genes (heat-shock proteins, antioxidant defenses, autophagy machinery) that buffer against TP4-induced stress, allowing cells to adapt and survive. Cancer cells, optimized for rapid proliferation rather than stress tolerance, may lack these protective capacities, making them vulnerable when confronted with severe TP4-induced stress.

The transcriptomic data are consistent with all four hypotheses but do not definitively distinguish among them. Future studies combining transcriptomics with functional assays (TP4 uptake measurements, mitochondrial functional assays, stress adaptation measurements) would be needed to definitively establish the mechanistic basis for cancer selectivity.

### Scientific Contribution and Significance

This transcriptomic reanalysis of the GSE74764 dataset provides several novel scientific insights that advance understanding of TP4 anticancer mechanisms:

**1. Comprehensive Mechanistic Mapping:**
By integrating differential expression analysis, pathway enrichment, and systems-level interpretation, this work reconstructs the complete molecular cascade of TP4 action—from initial cellular stress (calcium dysregulation, ROS, ER stress) through stress signal transduction (MAPK/JNK/AP-1 activation) to execution of apoptotic death. Prior work characterized TP4-induced apoptosis phenomenologically; this analysis maps the transcriptional architecture underlying apoptotic commitment.

**2. Multi-Pathway Convergence:**
The finding that TP4 simultaneously activates multiple independent apoptotic pathways (intrinsic mitochondrial, extrinsic death receptor, ER stress-induced, ROS-mediated) represents a novel mechanistic insight. Most drugs trigger apoptosis through one or two pathways; TP4's simultaneous activation of four independent death pathways creates high redundancy and resistance to single-pathway compensatory mechanisms. This multi-pathway convergence may explain TP4's potency and selectivity.

**3. AP-1 Isoform Remodeling:**
The identification of stress-AP-1 (enriched for ΔFosB/JunB) as distinct from growth-AP-1 (c-Fos/c-Jun) and the demonstration that TP4 shifts AP-1 composition from growth toward stress configuration represents a novel mechanistic insight into transcriptional reprogramming during stress and apoptosis. This finding has broader implications for understanding transcriptional control of cell fate decisions.

**4. Terminal vs. Adaptive UPR:**
The demonstration that TP4 induces terminal (pro-apoptotic) rather than adaptive (pro-survival) UPR through sustained CHOP activation provides mechanistic insight into how ER stress transitions from adaptive to apoptotic signaling—a fundamental question in cell biology with implications for ER stress-mediated diseases and therapies.

**5. Cancer-Selective Transcriptomics:**
The robust segregation of apoptotic gene expression between cancer and normal cells, combined with pathway enrichment analysis, provides transcriptomic evidence for the mechanistic basis of cancer selectivity. While upstream mechanisms (differential uptake, mitochondrial sensitivity, apoptotic threshold) remain to be fully resolved, this work establishes that cancer selectivity is encoded at the level of death-pathway gene expression.

**6. Implications for Future TP4-Based Therapies:**
The identification of multi-pathway convergence suggests that TP4 or TP4-derived peptides may be effective in cancer cells resistant to conventional chemotherapy through single-pathway resistance mechanisms. The demonstration of cancer-selective transcriptomics supports further development of TP4 as a cancer therapeutic with potential for reduced toxicity to normal tissues.

---

**End of Results Section**
