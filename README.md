# A Reproducibility-Focused Multi-Omics Latent-Factor Framework for Prognostic Modeling in Hepatocellular Carcinoma

**Author:** Abdullah I Ali  
**Affiliation:** Independent Researcher, Egypt  
**ORCID:** [0009-0005-6061-5562](https://orcid.org/0009-0005-6061-5562)  
**Corresponding Author Email:** abdallahhashem832@gmail.com  
**Additional Academic Contact:** abdullah416517@fsc.bu.edu.eg  

---

## Overview

This repository serves as the official **Model and Results Reproducibility Package** (Mode A) for the preprint study evaluating multi-modal prognostic modeling in Hepatocellular Carcinoma (TCGA-LIHC).

It distributes the sanitized, pre-trained multi-omics model bundle, exact feature specifications, aggregate cross-validation and holdout benchmark results, and production projection routines. To protect patient privacy, raw patient-level sequencing matrices and individual clinical survival records are not redistributed. Researchers wishing to perform end-to-end raw data reconstruction can acquire the source files directly from the NCI Genomic Data Commons (GDC) as documented in `docs/data_access.md`.

The study evaluates whether unsupervised multi-omics integration via Multi-Omics Factor Analysis v2 (MOFA2) can extract reproducible, prognostic latent representations across RNA sequencing, DNA methylation, copy number variation, and somatic mutations without data leakage.

> [!NOTE]
> **Study Reporting & Assessment Status:**
> - Reporting methodology evaluated against the TRIPOD (Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis) statement checklist; selected items are qualified according to study design (e.g., formal calibration curves were not estimated due to limited holdout death event counts).
> - All performance estimates reflect strictly out-of-fold cross-validation or projection into an untouched, locked internal holdout.
> - This work represents exploratory translational methodology and is **not approved, ready, or intended for clinical decision-making or bedside patient care**. Independent external cohort validation remains necessary.

---

## Study Design & Cohort Architecture

- **Initial Clinical Universe:** $N = 377$ patients diagnosed with primary hepatocellular carcinoma.
- **Pre-Analysis Exclusions ($N = 27$):** 5 non-primary tumor samples; 6 lacking follow-up time; 2 without primary pathology; 9 missing survival endpoints; 4 missing RNA/methylation assays; 1 extreme methylation outlier.
- **Primary Modeling Cohort ($N = 350$):**
  - **Development Cohort ($N = 285$):** 103 death events (36.14% event rate), 182 censored. Used for 5-repeat 5-fold nested cross-validation and final frozen model estimation.
  - **Locked Internal Holdout ($N = 65$):** 19 death events (29.23% event rate), 46 censored. Untouched throughout all model design, hyperparameter selection, and feature optimization.

---

## Multi-Omics Input Specifications

A total of **11,402 multi-modal features** were selected across four molecular layers:
1. **Gene Expression (RNA-seq):** 902 high-variance protein-coding genes (TMM-normalized, log2-CPM transformed, ComBat parametric batch adjusted for Tissue Source Site).
2. **DNA Methylation (Illumina 450K):** 5,000 top variable CpG probes (cross-reactive probes removed per Chen et al. 2013; ComBat batch adjusted).
3. **Copy Number Variation (Affymetrix SNP 6.0):** 5,000 high-variance genes derived from **Affymetrix SNP 6.0 precomputed segment data with length-weighted mean segment-to-gene aggregation**.
4. **Somatic Mutations (SNV):** 500 genes with highest mutation frequency filtered to nine functional classes: *Missense_Mutation, Nonsense_Mutation, Frame_Shift_Del, Frame_Shift_Ins, Splice_Site, In_Frame_Del, In_Frame_Ins, Nonstop_Mutation, Translation_Start_Site*.

---

## Methodological Innovations & Mathematical Formulation

### 1. Leakage-Free In-Fold Preprocessing
All batch adjustments, probe filtering, variance filtering, scaling, and imputation were conducted strictly inside each training fold to eliminate optimistic bias.

### 2. Regularized Out-of-Sample Latent Factor Projection
To project new samples into the frozen latent space without retraining MOFA2, a regularized joint ridge projection operator was derived:
$$\hat{z}_{\text{new}} = (W^T W + \lambda_{\text{ridge}} I)^{-1} W^T (x_{\text{new}} - \mu)$$
Setting $\lambda_{\text{ridge}} = 1.0$ regularizes the joint Gram matrix ($W^T W$), **providing empirical numerical stability in this study** (observed Gram condition number had a median of 773.4 and a maximum of 1632.9 across all 25 outer folds, with zero infinite or degenerate factor scores).

### 3. Penalized Survival Modeling
Factors were modeled using an ElasticNet-penalized Cox proportional hazards model (mixing parameter $\alpha = 0.5$; optimal shrinkage $\lambda_{\text{min}} = 0.08055$).

---

## Key Locked Scientific Results

- **Nested Cross-Validation Performance (25 runs):**
  - Mean $C$-index: **0.5875** (SD = 0.0707).
  - 5-Repeat Pooled $C$-index: **0.5844** (SD = 0.0164).
  - LASSO CV benchmark: 0.5910 (SD = 0.0678).
- **Locked Holdout Validation ($N = 65$):**
  - Harrell's $C$-index: **0.6197** (95% CI: 0.4785–0.7405, SE = 0.0573).
  - Time-Dependent Brier Scores: 1-Year = **0.0859**, 3-Year = **0.2180**, 5-Year = **0.2744**.
  - Integrated Brier Score (IBS): **0.1928** (null Kaplan-Meier reference: 0.2036).
- **Clinical Incremental Value:**
  - Clinical Stage alone: $C = 0.5974$.
  - Multi-Omics alone: $C = 0.5875$.
  - Combined Clinical + Omics: $C = 0.6120$ (Increment $\Delta C = +0.0146$, paired $t$-test $p = 0.2717$; not statistically conclusive).
- **Latent Factor 2 Biological Interpretation:**
  - Protective prognostic association: $\text{HR} = 0.7751$ (95% CI: 0.6770–0.8874, $p = 0.000224$).
  - High cross-fold recovery: Mean $|r| = 0.9784$, 100% selection frequency, 96% exact identification.
  - Enriched for differentiated mature hepatocyte metabolic pathways (xenobiotic metabolism, fatty acid degradation, bile acid synthesis).

---

## Repository Structure

```text
├── code/                   # Normalized reproduction scripts
│   ├── 01_preprocess_rna.R
│   ├── 02_preprocess_methylation.R
│   ├── 03_preprocess_cnv.R
│   ├── 04_preprocess_snv.R
│   ├── 05_generate_seed_fold_matrix.R
│   ├── 06_run_nested_cv.R
│   ├── 07_train_and_evaluate_pipeline.R
│   ├── 08_compute_pooled_cindex.R
│   ├── build_master_manuscript.py
│   └── verify_all_manuscript_facts.py
│   └── 09_apply_frozen_model.R
├── manifests/              # Patient identifiers (barcodes only)
│   ├── TRUE_DEV_285_IDS.csv
│   └── FINAL_TEST_IDS.csv
├── model/                  # Sanitized learned parameters
│   ├── FINAL_285_FROZEN_MODEL_BUNDLE_SANITIZED.rds
│   ├── MODEL_BUNDLE_MANIFEST.csv
│   ├── projection_operator.rds
│   └── elasticnet_coefficients.csv
├── data/references/        # Reference files and attrition ledgers
│   ├── Chen_2013_cross_reactive_probes.csv
│   ├── gene_coords_hg38.rds
│   ├── TRUE_DEV_285_SEED_FOLD_MATRIX.csv
│   └── GATE_2A_COHORT_ATTRITION.csv
├── results/                # Aggregate summary results tables
├── figures/                # Publication figures (main & supplementary)
│   ├── main/
│   └── supplementary/
└── docs/                   # Complete manuscript and evidence maps
    ├── MANUSCRIPT_MASTER.md
    ├── MANUSCRIPT_EVIDENCE_MAP.md
    ├── data_access.md
    ├── PUBLIC_RELEASE_MANIFEST.md
    ├── reproducibility.md
    └── manuscript_sections/
```

> [!IMPORTANT]
> **Data Privacy Notice on Patient-Level Predictions:**  
> In adherence to privacy-by-design principles, raw patient-level prediction logs (`nested_cv_285_patient_oof_predictions.csv`, `FINAL_285_TRUE_HOLDOUT_PREDICTIONS.csv`) are not redistributed in this public release. All aggregate fold metrics, summary tables, projection weights, and deterministic evaluation scripts are provided in full.

---

## Software & Environment Requirements

Verified computational stack:
- **R Version:** 4.4.1 (or 4.4+)
- **Python Version:** 3.14+ (with `mofapy2` 0.7.5 and `numpy`)
- **Key R Packages:** `MOFA2` (>= 1.16.0), `glmnet` (>= 4.1-10), `survival` (>= 3.8-3), `edgeR` (>= 4.4.2), `sva` (>= 3.54.0), `data.table` (>= 1.17.0).

---

## Citation & License

- **Code:** Licensed under the [MIT License](LICENSE) (Author confirmation pending).
- **Documentation, Results & Figures:** Licensed under [Creative Commons Attribution 4.0 International (CC-BY-4.0)](LICENSE).
- **Citation:** See [`CITATION.cff`](CITATION.cff) for bibtex and machine-readable metadata.
