# Reproducibility Framework & Execution Guidelines

This repository supports a multi-tiered reproducibility model for the manuscript:
**"A Reproducibility-Focused Multi-Omics Latent-Factor Framework for Prognostic Modeling in Hepatocellular Carcinoma"**

To uphold scientific transparency while strictly respecting data privacy principles (avoiding unauthorized redistribution of patient-level molecular matrices or uncurated clinical records), reproducibility is partitioned into three distinct, clearly defined levels.

---

## Reproducibility Hierarchy

```text
+-------------------------------------------------------------------------------+
| LEVEL 1: INSPECTION & VERIFICATION OF PUBLISHED AGGREGATE RESULTS             |
|   Inputs: Distributed results/ tables and figures/                            |
|   Goal: Audit reported metrics, variance tables, and factor statistics        |
|   Execution: Instantaneous; zero raw data or computational dependencies       |
+-------------------------------------------------------------------------------+
                                      |
                                      v
+-------------------------------------------------------------------------------+
| LEVEL 2: APPLICATION OF THE FROZEN MODEL TO COMPATIBLE MULTI-OMICS DATA       |
|   Inputs: model/FINAL_285_FROZEN_MODEL_BUNDLE_SANITIZED.rds                   |
|   Goal: Project new patient profiles into frozen MOFA space and predict risk   |
|   Execution: Pure R script (code/09_apply_frozen_model.R); lightweight        |
+-------------------------------------------------------------------------------+
                                      |
                                      v
+-------------------------------------------------------------------------------+
| LEVEL 3: FULL END-TO-END RECONSTRUCTION FROM GDC SOURCE ARCHIVES              |
|   Inputs: Raw GDC BAMs/counts/IDATs + GDC clinical tables                     |
|   Goal: Rerun in-fold preprocessing, ComBat, MOFA2 factorization, nested CV   |
|   Execution: code/01_preprocess_rna.R through code/07_train_pipeline.R        |
|   Requirements: External download of ~2.8 GB raw TCGA-LIHC data from GDC      |
+-------------------------------------------------------------------------------+
```

---

## Level 1: Verification of Published Aggregate Results

Researchers who wish to independently verify the numerical findings, statistical claims, factor recovery statistics, and figure panels can do so directly using the files distributed within this repository:

1. **Model Discrimination & Holdout Metrics:**
   - Inspect [`results/FINAL_285_TRUE_HOLDOUT_METRICS.csv`](../results/FINAL_285_TRUE_HOLDOUT_METRICS.csv) to verify locked holdout $C = 0.6197$, $\text{IBS} = 0.1928$, and time-dependent Brier scores ($1\text{Y}: 0.0859, 3\text{Y}: 0.2180, 5\text{Y}: 0.2744$).
2. **Cross-Validation Benchmarks (25 runs):**
   - Inspect [`results/nested_cv_285_fold_metrics.csv`](../results/nested_cv_285_fold_metrics.csv) for fold-by-fold metrics across all 5 repeats $\times$ 5 folds (ElasticNet mean $C = 0.5875$, LASSO mean $C = 0.5910$).
3. **Clinical Incremental Value:**
   - Inspect [`results/FINAL_CLINICAL_BENCHMARK_METRICS.csv`](../results/FINAL_CLINICAL_BENCHMARK_METRICS.csv) for Clinical ($C = 0.5974$), Multi-Omics ($C = 0.5875$), and Combined ($C = 0.6120, \Delta C = +0.0146, p = 0.2717$).
4. **Variance Decomposition:**
   - Inspect [`results/FINAL_MOFA_FACTOR_SUMMARY.csv`](../results/FINAL_MOFA_FACTOR_SUMMARY.csv) for modality-specific variance explained (CNV: 66.27%, Meth: 58.61%, RNA: 27.55%, SNV: 0.0892%).
5. **Manuscript Fact Audit:**
   - Run `python code/verify_all_manuscript_facts.py` to audit all 32 locked metrics and certify zero red-flag claim violations across the consolidated manuscript draft.

---

## Level 2: Application of Frozen Model to New Profiles

The pre-trained multi-omics prognostic model is distributed as a sanitized, learned parameter artifact:
`model/FINAL_285_FROZEN_MODEL_BUNDLE_SANITIZED.rds`

### Mathematical Formulation
New samples are projected into the frozen latent space using the regularized joint ridge projection operator ($\lambda_{\text{ridge}} = 1.0$):
$$\hat{z}_{\text{new}} = (W^T W + \lambda_{\text{ridge}} I)^{-1} W^T (x_{\text{new}} - \mu)$$

Prognostic risk scores are computed from the locked ElasticNet coefficients:
$$\hat{\eta} = \hat{z}_{\text{new}}^T \beta$$

### Execution
Run the standalone pure R engine:
```bash
Rscript code/09_apply_frozen_model.R
```
This script requires **zero Python dependencies** and operates purely through matrix operations in base R.

---

## Level 3: Full End-to-End Raw Reconstruction

> [!IMPORTANT]
> **Exact holdout metric recalculation ($C = 0.6197, \text{IBS} = 0.1928$) requires obtaining the raw TCGA-LIHC source files from GDC.**
> 
> Because raw molecular matrices are excluded from this repository, executing full upstream preprocessing and holdout scoring requires:
> 1. Downloading the raw TCGA-LIHC inputs using the instructions in [`docs/data_access.md`](data_access.md).
> 2. Partitioning the development ($N = 285$) and holdout ($N = 65$) sets using [`manifests/TRUE_DEV_285_IDS.csv`](../manifests/TRUE_DEV_285_IDS.csv) and [`manifests/FINAL_TEST_IDS.csv`](../manifests/FINAL_TEST_IDS.csv).
> 3. Running the sequential preprocessing pipelines:
>    - `code/01_preprocess_rna.R`
>    - `code/02_preprocess_methylation.R`
>    - `code/03_preprocess_cnv.R`
>    - `code/04_preprocess_snv.R`
> 4. Executing the master pipeline: `code/07_train_and_evaluate_pipeline.R`.
