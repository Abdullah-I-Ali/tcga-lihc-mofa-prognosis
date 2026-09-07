# SECTION 1: TITLE AND STRUCTURED ABSTRACT

## Title
**A Reproducibility-Focused Multi-Omics Latent-Factor Framework for Prognostic Modeling in Hepatocellular Carcinoma**

## Authors
[Author List / Affiliations Redacted for Preprint Submission]

---

## Abstract

### Background
High-dimensional genomic prognostic modeling is vulnerable to optimistic performance estimates when preprocessing, feature selection, or model development are not fully separated from evaluation. In hepatocellular carcinoma (HCC), multi-modal molecular profiling offers the potential to capture biological heterogeneity, yet the extent to which multi-omics latent representations add prognostic information beyond established clinical variables also remains uncertain.

### Methods
We developed and evaluated an integrative prognostic modeling framework using primary tumor profiles from The Cancer Genome Atlas Liver Hepatocellular Carcinoma (TCGA-LIHC) project. From an initial cohort of 377 patients, quality control and complete-case filtering across four modalities (transcriptomics, DNA methylation, copy number variation, and somatic mutations) yielded 350 eligible patients. A locked internal holdout of 65 patients (19 deaths, 46 censored) was quarantined before model development and remained untouched until final evaluation. The remaining 285 patients (103 deaths, 182 censored; median follow-up: 20.3 months) constituted the development cohort. Within the development cohort, a repeated 5-fold outer nested cross-validation protocol (5 repeats $\times$ 5 outer folds = 25 runs) was executed with strict in-fold encapsulation of all preprocessing steps: empirical Bayes batch correction with variance clamping ($\hat{\delta} \in [0.25, 4.0]$), modality-specific feature selection (902 RNA genes, 5,000 methylation CpGs, 5,000 CNV genes, and 500 somatic mutation events; 11,402 total features), and unsupervised Multi-Omics Factor Analysis (MOFA2; 15 latent factors). Out-of-sample factor scores were estimated using a joint ridge-regularized least-squares operator ($\lambda = 1.0$) on frozen training loadings rather than native MOFA2 prediction functions. Regularized Cox models (ElasticNet-Cox, $\alpha = 0.5$; LASSO-Cox, $\alpha = 1.0$), Random Survival Forests, and XGBoost-AFT models were tuned via inner 5-fold cross-validation. The final production ElasticNet-Cox model was fitted on the full development cohort ($N=285$) and evaluated exactly once on the quarantined locked internal holdout ($N=65$). A pre-specified clinical benchmark compared clinical covariates (age, gender, AJCC stage, grade), multi-omics factors, and combined models across all 25 outer folds.

### Results
In nested cross-validation on the development cohort, regularized Cox models demonstrated moderate discrimination: LASSO-Cox achieved a mean out-of-fold Harrell's C-index of 0.5910 (SD = 0.0678; mean Integrated Brier Score [IBS] = 0.2041), while ElasticNet-Cox achieved a mean C-index of 0.5875 (SD = 0.0707; pooled C-index = 0.5844; mean IBS = 0.2043; 1-year, 3-year, and 5-year Brier scores of 0.1422, 0.2275, and 0.2434, respectively). ElasticNet-Cox was selected as the primary production architecture to retain grouped correlated biological signals rather than for higher predictive discrimination. On the locked internal holdout (N=65), the model achieved a C-index of 0.6197 (95% bootstrap CI = 0.4785–0.7405; IBS = 0.1928; 1-year, 3-year, and 5-year Brier scores of 0.0859, 0.2180, and 0.2744, respectively). Dichotomizing holdout risk scores using the frozen development median threshold ($-0.0233$) yielded a hazard ratio of 0.834 (95% CI: 0.336–2.070, log-rank $p = 0.695$), with wide uncertainty around the dichotomized risk comparison. Biological audit of the 15 latent factors identified Factor 2 as the dominant prognostic latent axis (univariate $\text{HR} = 0.7751, 95\%\text{ CI: } 0.6770\text{--}0.8874, p = 0.000224, C = 0.6261$; ElasticNet $\beta = -0.1665$). Factor 2 proved highly reproducible across all 25 cross-validation runs (mean $|r| = 0.9784$, 100% selection frequency) and captured a transcriptomic contrast between mature hepatocyte differentiation/xenobiotic metabolism (positive loadings: *CYP3A4*, *CYP1A2*, *CYP8B1*; KEGG drug metabolism FDR = $1.93 \times 10^{-6}$) and mitotic cell-cycle proliferation (negative loadings: *MYBL2*, *SPP1*, *PEG10*; GO nuclear division FDR = $6.84 \times 10^{-25}$). In the clinical incremental-value benchmark, clinical covariates alone achieved a mean out-of-fold C-index of 0.5974 (SD = 0.0531), whereas multi-omics factors alone achieved a mean C-index of 0.5875 (SD = 0.0707). Combining clinical covariates with multi-omics factors yielded a mean C-index of 0.6120 (SD = 0.0676). The incremental gain of the combined model versus clinical staging alone was modest (mean $\Delta C = +0.0146$, nominal/exploratory paired $t$-test $p = 0.2717$; 5-repeat pooled $\Delta C = +0.0214$, bootstrap 95% CI: $-0.0372$ to $+0.0785$) and not statistically conclusive.

### Conclusions
Under rigorous cross-validation and out-of-sample projection, bulk multi-omics latent factors demonstrated moderate prognostic discrimination ($C \approx 0.59\text{--}0.62$). Multi-omics information appeared complementary to, rather than a replacement for, established clinical prognostic variables. The learned latent space captured an established differentiation-versus-proliferation biological axis with high stability across resamples, but the model is not intended for direct clinical decision-making. External validation in independent cohorts remains necessary to assess generalizability across diverse patient populations and assay platforms.

### Keywords
Hepatocellular carcinoma; Multi-omics integration; MOFA2; Prognostic modeling; Data leakage prevention; Internal validation; Survival analysis; Tumor differentiation.

---

## Compact Evidence Check: Section 1

| Metric / Fact Claimed in Abstract | Authoritative Source File | Source Value | Status in Draft |
| :--- | :--- | :--- | :---: |
| Development cohort $N=285$ (103 deaths, 182 censored) | `TRUE_DEV_285_IDS.csv`; `TRUE_DEV_285_EVENT_VECTOR.csv` | $N=285$, 103 deaths, 182 censored | **VERIFIED** |
| Holdout cohort $N=65$ (19 deaths, 46 censored) | `FINAL_TEST_IDS.csv`; `FINAL_285_TRUE_HOLDOUT_METRICS.csv` | $N=65$, 19 deaths, 46 censored | **VERIFIED** |
| Initial universe $N=377$, eligible complete cases $N=350$ | `GATE_2A7_FINAL_ACCOUNTING.csv` | $377 \to 350$ ($27$ excluded) | **VERIFIED** |
| Feature counts: 902 RNA, 5,000 Meth, 5,000 CNV, 500 SNV (11,402) | `GATE_2B5B_FINAL_MOFA_PILOT_AUDIT.md` | $902 + 5,000 + 5,000 + 500 = 11,402$ | **VERIFIED** |
| ComBat variance clamping $\hat{\delta} \in [0.25, 4.0]$ | `evaluate_285_production_and_holdout.R` | $\min = 0.25, \max = 4.0$ | **VERIFIED** |
| Ridge projection penalty $\lambda = 1.0$ | `evaluate_285_production_and_holdout.R` | $\lambda_{\text{ridge}} = 1.0$ | **VERIFIED** |
| CV LASSO mean $C = 0.5910$, IBS $= 0.2041$ | `FINAL_285_NESTED_CV_AUDIT.md` | Mean $C = 0.5910$, IBS $= 0.2041$ | **VERIFIED** |
| CV ElasticNet mean $C = 0.5875$, pooled $C = 0.5844$, IBS $= 0.2043$ | `FINAL_285_NESTED_CV_AUDIT.md` | Mean $C = 0.5875$, Pooled $0.5844$, IBS $0.2043$ | **VERIFIED** |
| ElasticNet not called "best-performing" | Abstract text | LASSO noted as 0.5910 vs 0.5875; ElasticNet chosen for stability | **VERIFIED** |
| Holdout Harrell's C-index $= 0.6197$ ($95\%\text{ CI: } 0.4785\text{--}0.7405$) | `FINAL_285_TRUE_HOLDOUT_METRICS.csv` | $0.6196769$ ($0.4784796\text{--}0.7404797$) | **VERIFIED** |
| Holdout IBS $= 0.1928$, Brier scores (0.0859, 0.2180, 0.2744) | `FINAL_285_TRUE_HOLDOUT_METRICS.csv` | IBS $= 0.19277$; 1y $0.0859$; 3y $0.2180$; 5y $0.2744$ | **VERIFIED** |
| Holdout KM binarization $\text{HR} = 0.834, p = 0.695$ | `FINAL_285_TRUE_HOLDOUT_METRICS.csv` | $\text{HR} = 0.83397$, log-rank $p = 0.69508$ | **VERIFIED** |
| Factor 2 univariate $\text{HR} = 0.7751, p = 0.000224, C = 0.6261$ | `FINAL_MOFA_FACTOR_SUMMARY.csv` | $\text{HR} = 0.7751, p = 0.0002236, C = 0.6261$ | **VERIFIED** |
| Factor 2 ElasticNet $\beta = -0.1665$ | `FINAL_MOFA_FACTOR_SUMMARY.csv` | $\beta = -0.166463$ | **VERIFIED** |
| Factor 2 CV stability $|r| = 0.9784$, 100% selection | `FINAL_FACTOR_STABILITY_AUDIT.csv` | Mean $|r| = 0.9784$, selection $= 100\%$ | **VERIFIED** |
| Factor 2 KEGG Drug metabolism FDR $= 1.93 \times 10^{-6}$ | `FINAL_PATHWAY_ENRICHMENT_AUDIT.csv` | $\text{FDR} = 1.930467 \times 10^{-6}$ | **VERIFIED** |
| Factor 2 GO Nuclear division FDR $= 6.84 \times 10^{-25}$ | `FINAL_PATHWAY_ENRICHMENT_AUDIT.csv` | $\text{FDR} = 6.837305 \times 10^{-25}$ | **VERIFIED** |
| Clinical benchmark: Clinical $C = 0.5974$, Omics $C = 0.5875$, Comb $C = 0.6120$ | `FINAL_CLINICAL_BENCHMARK_METRICS.csv` | Clin $0.5974$, Omics $0.5875$, Comb $0.6120$ | **VERIFIED** |
| Increment $\Delta C = +0.0146, p = 0.2717$; pooled $\Delta C = +0.0214$ | `FINAL_CLINICAL_INCREMENTAL_VALUE_AUDIT.md` | $\Delta C = +0.0146, p = 0.2717$; pooled $+0.0214$ | **VERIFIED** |
| Repeat 1 bootstrap 95% CI: $-0.0372$ to $+0.0785$ | `FINAL_CLINICAL_INCREMENTAL_VALUE_AUDIT.md` | $[-0.0372, +0.0785]$ | **VERIFIED** |

### Prohibited Red-Flag Word Check: Section 1
- `first` / `first-ever`: None
- `best` / `superior`: None
- `breakthrough` / `state-of-the-art`: None
- `clinically useful` / `clinically ready` / `ready for bedside`: None
- `external validation`: Prohibited use checked; text explicitly states "external validation in independent cohorts remains necessary".
- `invariant` / `causal` / `drives`: Prohibited causal language checked; Factor 2 termed "dominant prognostic latent axis", not causal driver.
- `0.6696`: None
- `IBS proves calibration`: None; framed as dynamic prediction error relative to null baselines.
- `ElasticNet best-performing`: None; explicitly acknowledged LASSO had higher mean C-index.

### Unresolved Statements in Section 1: None.
Section 1 is fully certified and internally consistent.
