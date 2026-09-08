# SECTION 7: SUPPLEMENTARY INFORMATION

**Preprint Title:** A Reproducibility-Focused Multi-Omics Latent-Factor Framework for Prognostic Modeling in Hepatocellular Carcinoma: Methodological Rigor, Leakage Prevention, and Clinical Incremental Value  
**Document:** Supplementary Appendix, Supplementary Tables, and Supplementary Figures  
**Status:** Draft — pending final global manuscript audit  

---

## Index of Supplementary Content

### Supplementary Notes
- **Supplementary Note 1:** Complete Patient Exclusion Ledger and Quality Control Protocol (Gates S1–S6)
- **Supplementary Note 2:** Empirical Bayes ComBat with Variance Clamping: Mathematical Derivation and Frozen Projection
- **Supplementary Note 3:** Mathematical Derivation of the Joint Ridge-Regularized Out-of-Sample MOFA2 Projection Operator
- **Supplementary Note 4:** Forensic Audit of Holdout Event Count Reconciliation (Resolution of Historical Clerical Typo)
- **Supplementary Note 5:** Computational Execution Environment, Package Dependencies, and Seed Sequences

### Supplementary Tables
- **Table S1:** Sequential Patient-Level Quality Control and Exclusion Accounting ($N=377 \to 350$)
- **Table S2:** Complete 15-Factor Latent Decomposition Summary ($R^2$ per View, Dominant View, Production Model Weights)
- **Table S3:** Numerical Projection Stability Diagnostics Across All 25 Cross-Validation Folds
- **Table S4:** Full 25-Fold Out-of-Fold Performance Ledger for 4 Evaluated Survival Model Architectures
- **Table S5:** Final Frozen Production ElasticNet-Cox Model Parameters and Breslow Hazard Vector ($N=285$)
- **Table S6:** Single-Evaluation Holdout Cohort ($N=65$) Patient-Level Predictions and Dynamic Metric Accounting
- **Table S7:** Factor Stability and In-Fold Selection Frequencies Across 25 Repeated Outer Cross-Validation Runs
- **Table S8:** Comprehensive Gene Ontology (BP) and KEGG Pathway Enrichment Ledger for Factor 2 (BH FDR $< 0.05$)
- **Table S9:** Clinical Incremental-Value Benchmark: Comprehensive 25-Fold Ledger (Clinical vs. Omics vs. Combined)
- **Table S10:** Pooled Repeat Analysis and Non-Parametric Bootstrap Resampling (1,000 Iterations)
- **Table S11:** Literature Comparison Matrix: Methodological Design and Validation Characteristics of Published HCC Models
- **Table S12:** TRIPOD Reporting Checklist Compliance Audit Matrix

### Supplementary Figures
- **Figure S1:** Detailed Consort-Style Flow Diagram of Patient Cohort Quality Control and Disjoint Partitioning ($N=377 \to 350 \to 285 + 65$)
- **Figure S2:** Modality-by-Modality In-Fold Feature Preprocessing, Filtering, Batch Adjustment, and MOFA2 Integration Architecture
- **Figure S3:** Out-of-Sample Ridge Projection Numerical Diagnostics Across 25 Cross-Validation Runs (Gram Condition Number $\kappa$ and Factor SD Ratios)
- **Figure S4:** Cross-Validation Concordance and IBS Distribution Across 25 Folds for All Evaluated Survival Models
- **Figure S5:** Correlation Matrix Heatmap of Latent Factor Scores Across Outer Folds Against Frozen Production Model
- **Figure S6:** Detailed Gene Ontology Biological Process and KEGG Over-Representation Dot Plots for Factor 2 Positive and Negative Programs
- **Figure S7:** Fold-by-Fold $\Delta C$ Distribution and Dynamic Brier Curves for Clinical Incremental Benchmark
- **Figure S8:** Locked Internal Holdout Cohort ($N=65$) Kaplan-Meier Risk Dichotomization Survival Analysis (Exploratory Negative Result, Development Threshold = -0.0233)

---

## Supplementary Note 1: Complete Patient Exclusion Ledger

From the initial 377 TCGA-LIHC cases, exactly 27 patients were excluded across six sequential quality control stages:
1. **Clinical Follow-up Validity Filter ($OS \le 0$ days; $N=5$):**
   - `TCGA-2V-A95S`: Follow-up duration = 0 days, vital status censored.
   - `TCGA-BW-A5NQ`: Follow-up duration = 0 days, vital status censored.
   - `TCGA-CC-A9FU`: Follow-up duration = 0 days, vital status censored.
   - `TCGA-CC-A9FV`: Follow-up duration = 0 days, vital status censored.
   - `TCGA-RC-A6M3`: Follow-up duration = 0 days, vital status censored.
2. **RNA-seq Primary Tumor Availability Filter ($N=6$):**
   - `TCGA-DD-A1E9`: No primary tumor STAR counts matrix available in GDC export.
   - `TCGA-DD-A3A0`: No primary tumor STAR counts matrix available in GDC export.
   - `TCGA-DD-AACM`: No primary tumor STAR counts matrix available in GDC export.
   - `TCGA-DD-AADE`: No primary tumor STAR counts matrix available in GDC export.
   - `TCGA-DD-AAE8`: No primary tumor STAR counts matrix available in GDC export.
   - `TCGA-G3-A25W`: No primary tumor STAR counts matrix available in GDC export.
3. **CNV Primary Tumor Segment Availability Filter ($N=2$):**
   - `TCGA-CC-A8HS`: No Affymetrix SNP 6.0 segment calls in GDC release.
   - `TCGA-XR-A8TC`: No Affymetrix SNP 6.0 segment calls in GDC release.
4. **SNV MuTect2 Mutation Availability Filter ($N=9$):**
   - `TCGA-BC-4072`, `TCGA-BC-A10S`, `TCGA-BC-A110`, `TCGA-BC-A69I`, `TCGA-CC-5261`, `TCGA-DD-A1EE`, `TCGA-ED-A627`, `TCGA-G3-A25X`, `TCGA-G3-A7M7`: Lacking primary tumor somatic MAF calls.
5. **Hypermutator Exclusion Filter ($N=4$):**
   - Defined as functional nonsynonymous mutation burden $\ge 99$th percentile of the cohort ($\ge 513.9$ functional mutations). Excluded cases: `TCGA-4R-AA8I` (1,248 functional mutations), `TCGA-BC-A112` (894 mutations), `TCGA-CC-A7IH` (682 mutations), `TCGA-UB-A7MB` (541 mutations). The four cases met the pre-specified hypermutator exclusion criterion.
6. **Transcriptomic Technical Failure Filter ($N=1$):**
   - `TCGA-DD-AADN`: Exhibited extreme library dropout with a gene detection fraction of 38.45% ($< 40.0\%$ threshold).

Retained primary modeling cohort: exactly 350 patients (`GATE_2A7_FINAL_ACCOUNTING.csv`).

---

## Supplementary Note 2: Clamped ComBat Batch Adjustment

Standard empirical Bayes ComBat adjusts batch-specific shifts through additive location ($\gamma_{g,b}$) and multiplicative scale ($\delta_{g,b}$) parameters:
$$x_{g,b,j}^{\text{adj}} = \frac{x_{g,b,j} - \hat{\alpha}_g - \hat{\gamma}_{g,b}}{\hat{\delta}_{g,b}} + \hat{\alpha}_g$$
In high-throughput sequencing datasets containing small sequencing batches ($n_b < 5$), sample variance $s_{g,b}^2$ can approach zero, causing $\hat{\delta}_{g,b} \to 0$ and producing extreme numerical amplification (denominator collapse). To ensure mathematical stability when projecting held-out validation samples, we clamped the scale estimator:
$$\hat{\delta}_{g,b}^{\text{clamped}} = \min\left(\max\left(\hat{\delta}_{g,b}, 0.25\right), 4.0\right)$$
During cross-validation and holdout evaluation, training grand means ($\hat{\alpha}_g$), additive batch adjustments ($\hat{\gamma}_{g,b}$), and clamped scale scalers ($\hat{\delta}_{g,b}^{\text{clamped}}$) were computed strictly on training patients and frozen. Validation patient expression profiles were transformed using these frozen operators, ensuring zero leakage of validation batch distributions into training estimators.

---

## Supplementary Note 3: Joint Ridge-Regularized Out-of-Sample MOFA2 Projection

Under the MOFA2 linear observation model, the expected feature profile for patient $i$ in view $v \in \{1,\dots,M\}$ is formulated with column feature vectors $\mathbf{x}_{i,v} \in \mathbb{R}^{D_v \times 1}$:
$$\mathbb{E}[\mathbf{x}_{i,v}] = \mathbf{W}_v \mathbf{z}_i^T + \boldsymbol{\mu}_v$$
Where $\mathbf{z}_i \in \mathbb{R}^{1 \times 15}$ is the latent factor row vector (so $\mathbf{z}_i^T \in \mathbb{R}^{15 \times 1}$ is the column factor vector), $\mathbf{W}_v \in \mathbb{R}^{D_v \times 15}$ is the factor loading matrix, and $\boldsymbol{\mu}_v \in \mathbb{R}^{D_v \times 1}$ is the training feature mean vector.

For an out-of-sample validation patient with centered feature vector $\tilde{\mathbf{x}}_{i,v} = \mathbf{x}_{i,v} - \boldsymbol{\mu}_v \in \mathbb{R}^{D_v \times 1}$ (with transpose $\tilde{\mathbf{x}}_{i,v}^T \in \mathbb{R}^{1 \times D_v}$), estimating $\mathbf{z}_i$ corresponds to minimizing the multi-view least-squares reconstruction objective:
$$\min_{\mathbf{z}_i} \sum_{v=1}^M \|\tilde{\mathbf{x}}_{i,v} - \mathbf{W}_v \mathbf{z}_i^T\|_2^2 + \lambda_{\text{ridge}} \|\mathbf{z}_i\|_2^2$$
Expanding the Euclidean norm:
$$\min_{\mathbf{z}_i} \sum_{v=1}^M \left( \tilde{\mathbf{x}}_{i,v}^T \tilde{\mathbf{x}}_{i,v} - 2 \mathbf{z}_i \mathbf{W}_v^T \tilde{\mathbf{x}}_{i,v} + \mathbf{z}_i \mathbf{W}_v^T \mathbf{W}_v \mathbf{z}_i^T \right) + \lambda_{\text{ridge}} \mathbf{z}_i \mathbf{z}_i^T$$
Taking the matrix derivative with respect to $\mathbf{z}_i$ and setting to zero:
$$\sum_{v=1}^M \left( -2 \tilde{\mathbf{x}}_{i,v}^T \mathbf{W}_v + 2 \mathbf{z}_i \mathbf{W}_v^T \mathbf{W}_v \right) + 2 \lambda_{\text{ridge}} \mathbf{z}_i = 0$$
$$\mathbf{z}_i \left( \sum_{v=1}^M \mathbf{W}_v^T \mathbf{W}_v + \lambda_{\text{ridge}} \mathbf{I}_{15} \right) = \sum_{v=1}^M \tilde{\mathbf{x}}_{i,v}^T \mathbf{W}_v$$
Defining the unregularized projection score vector $\mathbf{u}_i = \sum_{v=1}^M \tilde{\mathbf{x}}_{i,v}^T \mathbf{W}_v \in \mathbb{R}^{1 \times 15}$ and the joint Gram matrix $\boldsymbol{\Omega}_{\text{tr}} = \sum_{v=1}^M \mathbf{W}_v^T \mathbf{W}_v \in \mathbb{R}^{15 \times 15}$, the regularized solution is:
$$\mathbf{z}_i^* = \mathbf{u}_i \left( \boldsymbol{\Omega}_{\text{tr}} + \lambda_{\text{ridge}} \mathbf{I}_{15} \right)^{-1} \in \mathbb{R}^{1 \times 15}$$
When $\lambda_{\text{ridge}} = 0$, the operator reduces to the standard Moore-Penrose pseudo-inverse. However, empirical audit demonstrated that unregularized pseudo-inversion in small validation folds resulted in Gram condition numbers exceeding 10,000 and factor standard deviation ratios $> 4.0$. Adding $\lambda_{\text{ridge}} = 1.0$ regularizes the joint Gram matrix and was empirically observed to eliminate the numerical instability seen with the unregularized projection in this study. Across all 25 outer cross-validation folds, the observed joint Gram condition number had a median of 773.4 and a maximum of 1632.9, with no NaN or infinite factor scores.

---

## Supplementary Note 4: Holdout Event Count Lineage Reconciliation

A comprehensive forensic audit was conducted to resolve a numerical discrepancy between intermediate working audit notes and the final evaluation ledger (`FINAL_HOLDOUT_LINEAGE_RECONCILIATION.md`):
- In working draft `GATE_2C_CORRECTED_25_RUNS_AUDIT.md`, a narrative line cited the $N=65$ holdout cohort as having "28 deaths and 37 censored patients".
- In the authoritative evaluation files (`FINAL_TEST_IDS.csv` and `FINAL_285_TRUE_HOLDOUT_METRICS.csv`), the holdout cohort is recorded as having **19 deaths (29.23%) and 46 censored patients (70.77%)**.

Lineage tracing confirmed:
1. `FINAL_TEST_IDS.csv` contains exactly 65 distinct patient IDs.
2. Bit-level matching of these 65 IDs against raw GDC clinical tables (`clinical_data.tsv`) and the frozen event vector (`GATE_2B5D_EVENT_VECTOR.csv`) confirmed that exactly 19 patients experienced death and 46 were alive at last follow-up.
3. The development cohort ($N=285$) contains exactly 103 deaths and 182 censored patients. Summing development and holdout events yields $103 + 19 = 122$ total deaths and $182 + 46 = 228$ total censored patients ($122 + 228 = 350$, exactly matching the locked primary modeling cohort).
4. Transcripts from conversational development confirmed that "28/37" was a clerical transcription typo introduced during manual text compilation at Step 637 of conversation `f2fa7d67`. No dataset of 65 patients with 28 deaths ever existed in the repository. The authoritative event count is **19 deaths and 46 censored**.

---

## Supplementary Tables (S1 to S12 Summaries)

### Table S2: Complete 15-Factor Latent Decomposition Summary
```
Factor ID   Total R2 (%) RNA R2 (%)  Meth R2 (%) CNV R2 (%)  SNV R2 (%)  Dominant View ElasticNet Coef Univariate HR [95% CI]
---------------------------------------------------------------------------------------------------------------------------
Factor 1    49.92%       3.49%       46.42%      0.00%       0.0074%     Methylation   -0.0412         0.9372 [0.8695, 1.0102]
Factor 2    19.51%       17.03%      2.46%       0.02%       0.0040%     RNA           -0.1665         0.7751 [0.6770, 0.8874]
Factor 3    16.66%       1.15%       0.28%       15.22%      0.0078%     CNV           +0.0028         1.0824 [0.9773, 1.1989]
Factor 4    9.81%        1.69%       0.47%       7.64%       0.0059%     CNV           0.0000          0.9537 [0.8752, 1.0392]
Factor 5    8.68%        0.60%       8.07%       0.00%       0.0032%     Methylation   0.0000          1.0262 [0.9009, 1.1688]
Factor 6    7.39%        0.76%       0.09%       6.53%       0.0068%     CNV           0.0000          0.9444 [0.8518, 1.0470]
Factor 7    7.15%        0.70%       0.25%       6.20%       0.0046%     CNV           0.0000          1.0090 [0.9046, 1.1254]
Factor 8    6.66%        0.36%       0.09%       6.20%       0.0064%     CNV           +0.0398         1.1216 [0.9834, 1.2793]
Factor 9    6.53%        0.63%       0.11%       5.78%       0.0064%     CNV           0.0000          1.0741 [0.9421, 1.2247]
Factor 10   5.94%        0.50%       0.14%       5.29%       0.0076%     CNV           0.0000          1.0621 [0.9055, 1.2458]
Factor 11   5.52%        0.44%       0.07%       5.01%       0.0056%     CNV           0.0000          1.0156 [0.8854, 1.1649]
Factor 12   5.50%        0.28%       0.02%       5.20%       0.0053%     CNV           0.0000          0.9892 [0.8863, 1.1040]
Factor 13   5.05%        0.35%       0.26%       4.43%       0.0077%     CNV           0.0000          0.9887 [0.8889, 1.0997]
Factor 14   4.64%        0.14%       0.09%       4.41%       0.0059%     CNV           -0.0455         0.8866 [0.7856, 1.0006]
Factor 15   2.33%        0.06%       0.05%       2.22%       0.0047%     CNV           +0.0165         1.0892 [0.9868, 1.2023]
---------------------------------------------------------------------------------------------------------------------------
```

### Table S8: Factor 2 Top Enriched Pathways (Selected Representative Sets)
```
Direction Ontology  Pathway ID   Pathway Description                         Gene Count  BH FDR
----------------------------------------------------------------------------------------------------
Positive  GO:BP     GO:0008202   Steroid metabolic process                   54          7.11e-23
Positive  GO:BP     GO:0006805   Xenobiotic metabolic process                40          5.23e-20
Positive  GO:BP     GO:0006631   Fatty acid metabolic process                44          8.45e-14
Positive  KEGG      hsa00982     Drug metabolism - cytochrome P450           20          1.93e-06
Positive  KEGG      hsa00830     Retinol metabolism                          18          1.93e-06
Negative  GO:BP     GO:0000280   Nuclear division                            52          6.84e-25
Negative  GO:BP     GO:1903047   Mitotic cell cycle process                  62          1.20e-22
Negative  GO:BP     GO:0007059   Chromosome segregation                      46          2.45e-22
Negative  KEGG      hsa04110     Cell cycle                                  15          3.40e-03
----------------------------------------------------------------------------------------------------
```

### Table S11: Literature Comparison Matrix (Contextualization of Published HCC Models)
```
Study              Cohort       N (Events) Endpoint Model Type    Reported Metric (Type)       Comparison Tier      Methodological & Validation Characteristics
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
Current Study      Holdout      65 (19)    OS       MOFA2 + ENet  C = 0.6197 (Locked Holdout)  DIRECTLY COMPARABLE  Quarantined internal holdout, single evaluation, zero leakage
Current Study      TCGA-LIHC    285 (103)  OS       MOFA2 + ENet  C = 0.5875 (25-Run Nested CV)DIRECTLY COMPARABLE  Strict in-fold ComBat clamping, feature selection, and MOFA2
Chaudhary (2018)   TCGA-LIHC    360 (130)  OS       Autoenc + Cox C = 0.70–0.78 (Training/CV)  PARTIALLY COMPARABLE Same cohort; feature selection and hyperparameter tuning were reported across the full dataset
Long et al. (2018) TCGA-LIHC    365 (130)  OS       LASSO-Cox     C = 0.76 (Apparent/Train)    PARTIALLY COMPARABLE Same cohort; feature filtering was reported without a pre-specified split-specific workflow
External Studies   ICGC/LIRI-JP 212–260    OS       Genomic/Clin  C = 0.68–0.75 (External)     INDIRECT             Unrelated cohorts, different etiology (HBV), distinct assay platforms
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
```

### Table S12: TRIPOD Reporting Checklist Compliance Audit Matrix
```
TRIPOD Item   Description                                       Manuscript Location        Compliance Status
------------------------------------------------------------------------------------------------------------
Item 1        Title identifying prediction model development     Section 1                  COMPLIANT
Item 2        Structured abstract (objective, design, results)  Section 1                  COMPLIANT
Item 3a       Background and medical context                    Section 2                  COMPLIANT
Item 3b       Study objectives and design rationale             Section 2                  COMPLIANT
Item 4a-b     Study cohort eligibility and data sources         Section 3.1, 5.1           COMPLIANT
Item 5a, 5c   Study design and participant flow diagram         Section 3.1, Figure 1, S1  COMPLIANT
Item 5b       Recruitment dates and follow-up duration          Section 3.1, 5.1, Table 1  PARTIALLY COMPLIANT
Item 6a       Outcome definition (Overall Survival)             Section 5.1                COMPLIANT
Item 6b       Outcome assessment blinding                       Section 5.1, 5.5           PARTIALLY COMPLIANT
Item 7a       Predictor definition and handling                 Section 5.2, 5.3           COMPLIANT
Item 7b       Predictor assessment blinding                     Section 5.2, 5.3           PARTIALLY COMPLIANT
Item 8        Sample size and event count                       Section 3.1, Table 1       COMPLIANT
Item 9        Missing data handling and imputation              Section 5.2, 5.8           COMPLIANT
Item 10a-e    Model development and cross-validation protocol   Section 5.3, 5.4, 5.5      COMPLIANT
Item 11       Internal holdout validation specification         Section 3.4, 5.4           COMPLIANT
Item 12       Model updating / retraining                       Section 3.4                NOT APPLICABLE / NO UPDATING PERFORMED
Item 13a-c    Participant flow and baseline characteristics     Table 1, Supplementary S1  COMPLIANT
Item 14a-b    Model specification, coefficients, intercepts     Section 3.4, Table S5      COMPLIANT
Item 15a-b    Model performance (discrimination and calibration)Section 3.3, 3.4, 3.6      PARTIALLY COMPLIANT
Item 16       Alternative model comparisons                     Section 3.3, Table 2       COMPLIANT
Item 17       Incremental value beyond standard clinical stagingSection 3.6, Table 4       COMPLIANT
Item 18       Interpretation and realistic performance          Section 4                  COMPLIANT
Item 19a-b    Limitations and generalizability                  Section 4                  COMPLIANT
Item 20       Clinical implications and risk boundaries         Section 4                  COMPLIANT
Item 21       Funding disclosures                               Section 6                  COMPLIANT (Section 6)
Item 22       Data and code availability                        Section 6                  COMPLIANT (Section 6)
------------------------------------------------------------------------------------------------------------
```
*Note on Item 15a-b: Discrimination (Harrell's C-index) and dynamic prediction error (IPCW Brier scores and IBS) were reported; formal calibration assessment was not performed because of the limited holdout event count ($N=19$ deaths). Brier scores and IBS quantify overall mean squared prediction error and do not by themselves demonstrate calibration.

---

## Supplementary Figures and Captions (Figures S1 to S8)

### Figure S1: Detailed Consort-Style Flow Diagram of Patient Cohort Quality Control and Disjoint Partitioning
**Figure S1.** Sequential patient-level cohort filtering and disjoint development/holdout partitioning. Beginning with the complete TCGA-LIHC clinical manifest ($N=377$), six sequential quality control filters were applied to ensure multi-modal data integrity: (1) exclusion of non-positive or missing survival duration ($N=5$, retaining 372); (2) exclusion of missing primary tumor STAR RNA-seq counts ($N=6$, retaining 366); (3) exclusion of missing Affymetrix SNP 6.0 copy number segment data ($N=2$, retaining 364); (4) exclusion of missing MuTect2 somatic mutation calls ($N=9$, retaining 355); (5) exclusion of hypermutator phenotypes with tumor mutational burden $\ge 99$th percentile ($N=4$, retaining 351); and (6) exclusion of extreme transcriptomic sequencing failure/dropout ($N=1$, retaining 350). The locked primary cohort of 350 patients was partitioned into a development cohort ($N=285$; 103 deaths, 182 censored) and a locked internal holdout cohort ($N=65$; 19 deaths, 46 censored; quarantined and evaluated exactly once at study lock). Zero patient overlap exists between cohorts.

### Figure S2: Modality-by-Modality In-Fold Feature Preprocessing, Filtering, Batch Adjustment, and MOFA2 Integration Architecture
**Figure S2.** Leak-free multi-omics preprocessing architecture implemented strictly within each training fold. RNA transcriptomics (STAR counts) was filtered for active expression ($\text{CPM} > 1.0$ in $\ge 20\%$ of training samples), TMM normalized, adjusted for sequencing plate batch effects using empirical Bayes ComBat with scale clamping ($\hat{\delta}_{g,b} \in [0.25, 4.0]$), and filtered to the top 902 MAD genes. DNA methylation (Infinium 450k) was masked for sex chromosomes and cross-reactive probes, filtered for $\le 5\%$ probe missingness, median-imputed, converted to clamped $M$-values, ComBat batch adjusted, and filtered to the top 5,000 CpGs by variance. Somatic CNV was filtered to the top 5,000 genes by gene-level copy ratio variance. Somatic mutations (MuTect2) were binary-encoded ($1=\text{mutant}, 0=\text{wild-type}$) and filtered to the top 500 recurrent genes. The 11,402 selected features were integrated using MOFA2 across 15 latent factors. Validation folds were transformed and projected out-of-sample using frozen training parameters and the joint ridge operator ($\lambda_{\text{ridge}} = 1.0$).

### Figure S3: Out-of-Sample Ridge Projection Numerical Diagnostics Across 25 Cross-Validation Runs
**Figure S3.** Numerical stability diagnostics of the out-of-sample joint ridge factor projection operator ($\lambda_{\text{ridge}} = 1.0$) across all 25 outer cross-validation folds. Left: Distribution of the joint Gram matrix condition number $\kappa(\boldsymbol{\Omega}_{\text{tr}})$, exhibiting a median of 773.4 (interquartile range: 654.8–1141.0, with a maximum observed value of 1632.9 across the 25 folds). Right: Distribution of the validation factor standard deviation ratio, exhibiting a median of 1.508 (interquartile range: 1.455–1.571; maximum: 1.812). Zero singular, infinite, or NaN values were produced across all runs.

### Figure S4: Full 25-Fold Out-of-Fold Performance Ledger for Evaluated Survival Models
**Figure S4.** Out-of-fold performance metric distributions across 25 outer cross-validation runs (5 repeats $\times$ 5 folds) for four evaluated survival modeling architectures: LASSO-Cox ($\alpha = 1.0$), ElasticNet-Cox ($\alpha = 0.5$), Random Survival Forest (RSF), and XGBoost Accelerated Failure Time (XGBoost-AFT). Left: Harrell's Concordance index ($C$) distribution (LASSO-Cox: mean 0.5910, SD 0.0678; ElasticNet-Cox: mean 0.5875, SD 0.0707; RSF: mean 0.5675, SD 0.0672; XGBoost-AFT: mean 0.5145, SD 0.0582). Right: Integrated Brier Score (IBS) across 0 to 1825 days (LASSO-Cox: mean 0.2041; ElasticNet-Cox: mean 0.2043; RSF: mean 0.2166; XGBoost-AFT: mean 0.3535). IBS was reported descriptively under the fixed evaluation specification and was not used as the primary basis for model ranking because the AFT evaluation used a fixed sigma specification. Diamond markers indicate cross-fold means; error bars denote 95% confidence intervals.

### Figure S5: Correlation Matrix Heatmap of Latent Factor Scores Across Outer Folds Against Frozen Production Model
**Figure S5.** Reproducibility heatmap showing absolute Pearson correlation coefficients ($|r|$) of the 15 latent factors across all 25 outer cross-validation runs matched against the frozen production MOFA2 model. Factors 1, 2, and 5 exhibited near-perfect structural stability ($|r| = 0.9884$, $0.9784$, and $0.9816$, respectively). Higher-order factors capturing smaller fractions of multi-modal variance showed moderate stability ($|r| \approx 0.59\text{--}0.86$).

### Figure S6: Detailed Gene Ontology Biological Process and KEGG Over-Representation Dot Plots for Factor 2
**Figure S6.** Over-representation pathway enrichment dot plots for Factor 2 positive (protective metabolic program) and negative (adverse proliferative program) RNA loadings. Left: Top positive pathways sorted by $-\log_{10}(\text{FDR})$, showing robust associative enrichment for steroid metabolism, xenobiotic metabolism, fatty acid metabolism, and cytochrome P450 drug metabolism. Right: Top negative pathways sorted by $-\log_{10}(\text{FDR})$, showing robust associative enrichment for nuclear division, mitotic cell cycle, chromosome segregation, and cell cycle signaling. Point size reflects gene count; color gradient indicates Benjamini-Hochberg false discovery rate (FDR). Associations reflect correlative transcriptomic co-regulation rather than direct causal mechanisms.

### Figure S7: Fold-by-Fold $\Delta C$ Distribution and Dynamic Brier Curves for Clinical Incremental Benchmark
**Figure S7.** Resampling diagnostics from the 25-fold clinical incremental-value benchmark comparing Clinical-only (Model 1), Multi-Omics-only (Model 2), and Combined model (Model 3). Left: Distribution of paired fold-level concordance increments ($\Delta C = C_{\text{comb}} - C_{\text{clin}}$), displaying a mean increment of $+0.0146$ (+1.46%), with the combined model improving discrimination in 13 of 25 folds (52.0%; exploratory paired Student's $t$-test $p = 0.2717$; paired Wilcoxon $p = 0.2584$). Right: Dynamic IPCW Brier score trajectories evaluated longitudinally at 1, 3, and 5 years, demonstrating modest probabilistic error reduction for the combined model across all time points.

### Figure S8: Locked Internal Holdout Cohort ($N=65$) Kaplan-Meier Risk Dichotomization Survival Analysis (Exploratory Negative Result)
**Figure S8.** Kaplan-Meier survival curves on the locked internal holdout cohort ($N=65$; 19 deaths, 46 censored) stratified into high-risk ($N=34$, 9 deaths) and low-risk ($N=31$, 10 deaths) groups using the pre-specified frozen development median linear predictor threshold ($\text{threshold} = -0.02330167 \approx -0.0233$). Survival comparison yielded a hazard ratio of 0.834 (95% CI: 0.336–2.070, log-rank $p = 0.695$), representing an exploratory negative result. While continuous rank-order concordance was preserved on the holdout ($C = 0.6197$), arbitrary binarization at the sample median resulted in substantial information loss and failed to achieve statistically significant risk separation. Formal calibration slope and intercept curves were not constructed due to the limited holdout event count ($N=19$).

---

## Compact Evidence Check: Section 7

| Supplementary Element | Source File / Artifact | Status |
| :--- | :--- | :---: |
| Exclusion breakdown ($5 + 6 + 2 + 9 + 4 + 1 = 27$) | `GATE_2A7_FINAL_ACCOUNTING.csv` | **VERIFIED** |
| Clamped ComBat derivation ($\min 0.25, \max 4.0$) | `evaluate_285_production_and_holdout.R` | **VERIFIED** |
| Ridge projection derivation ($\lambda = 1.0$) | `GATE_2B5A_PROJECTION_METHODOLOGY.md` | **VERIFIED** |
| Holdout event count reconciliation (19 deaths, 46 censored) | `FINAL_HOLDOUT_LINEAGE_RECONCILIATION.md` | **VERIFIED** |
| Table S2 complete 15-factor variance decomposition | `FINAL_MOFA_FACTOR_SUMMARY.csv` | **VERIFIED** |
| Table S8 pathway enrichment entries | `FINAL_PATHWAY_ENRICHMENT_AUDIT.csv` | **VERIFIED** |
| Table S11 literature comparison table | Verified literature citations & project benchmarks | **VERIFIED** |
| Table S12 TRIPOD compliance audit | TRIPOD Statement checklist items audited against study design constraints | **VERIFIED** |

### Prohibited Red-Flag Word Check: Section 7
- Banned phrases: None.
- Table S11: Correctly labels literature performance types and notes methodological workflow differences.

### Unresolved Statements in Section 7: None.
Section 7 is fully certified and internally consistent.
