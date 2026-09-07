# Public Release Artifact Manifest

This document provides a comprehensive audit trail of all files prepared in the public staging repository (`tcga-lihc-mofa-prognosis-public`). Every staged artifact has been verified to ensure full scientific reproducibility while strictly excluding raw patient-level molecular matrices and individual-level clinical covariates.

---

## Summary of Staged Files by Category

| Relative Path | Category | Purpose | Public Release Status | Contains Patient Data? | Source Path in Analysis Workspace | Caveats / Notes |
| :--- | :--- | :--- | :---: | :---: | :--- | :--- |
| `README.md` | Documentation | Repository overview and reproducibility guide | **APPROVED** | **NO** | New Draft | Adheres to non-claim restrictions; specifies actual CNV aggregation |
| `LICENSE` | Legal | Dual MIT / CC-BY-4.0 licensing terms | **APPROVED** | **NO** | New Draft | Subject to author final confirmation |
| `CITATION.cff` | Metadata | Machine-readable citation record | **APPROVED** | **NO** | New Draft | Formatted for Abdullah I Ali |
| `.gitignore` | Config | Prevents accidental commit of raw matrices/logs | **APPROVED** | **NO** | New Draft | Explicitly blocks patient-level data |
| `.gitattributes` | Config | Git LFS tracking for `.rds`, `.pdf`, `.png` | **APPROVED** | **NO** | New Draft | Preserves lightweight Git history |
| `manifests/TRUE_DEV_285_IDS.csv` | Manifest | Development cohort patient barcodes ($N=285$) | **APPROVED** | **NO** (GDC barcodes only) | `TRUE_DEV_285_IDS.csv` | Single column `patient_id` |
| `manifests/FINAL_TEST_IDS.csv` | Manifest | Locked holdout patient barcodes ($N=65$) | **APPROVED** | **NO** (GDC barcodes only) | `FINAL_TEST_IDS.csv` | Sanitized: demographic columns (`age`, `gender`, `stage`) removed |
| `model/FINAL_285_FROZEN_MODEL_BUNDLE_SANITIZED.rds` | Model | Complete learned model and projection parameters | **APPROVED** | **NO** | Sanitized from `FINAL_285_FROZEN_MODEL_BUNDLE.rds` | Raw molecular data matrices (`@data`), sample barcodes, and factor scores stripped |
| `model/MODEL_BUNDLE_MANIFEST.csv` | Model | Machine-readable audit of bundle parameters | **APPROVED** | **NO** | Generated during sanitization | Documents dimensions and classes of retained objects |
| `model/projection_operator.rds` | Model | Decoupled linear algebra Gram projection operator | **APPROVED** | **NO** | `FINAL_FROZEN_MODEL_BUNDLE/projection_operator.rds` | Pure linear algebra: $(W^T W + \lambda I)^{-1}$ |
| `model/elasticnet_coefficients.csv` | Model | Locked factor regression coefficients | **APPROVED** | **NO** | `FINAL_FROZEN_MODEL_BUNDLE/elasticnet_coefficients.csv` | 6 non-zero Cox hazard weights |
| `data/references/Chen_2013_cross_reactive_probes.csv` | Reference | Quality control cross-reactive 450K probe list | **APPROVED** | **NO** | `Chen_2013_cross_reactive_probes.csv` | Literature probe list (Chen et al. 2013) |
| `data/references/gene_coords_hg38.rds` | Reference | Genomic coordinates for segment-to-gene mapping | **APPROVED** | **NO** | `gene_coords_hg38.rds` | hg38 gene annotation coordinates |
| `data/references/TRUE_DEV_285_SEED_FOLD_MATRIX.csv` | Reference | Exact 5-repeat 5-fold CV assignment matrix | **APPROVED** | **NO** | `TRUE_DEV_285_SEED_FOLD_MATRIX.csv` | Fold assignments (1 to 5) across 5 repeats |
| `data/references/GATE_2A_COHORT_ATTRITION.csv` | Reference | Step-by-step cohort exclusion ledger ($377 	o 350$) | **APPROVED** | **NO** | `GATE_2A_COHORT_ATTRITION.csv` | Aggregate counts and exclusion codes |
| `results/FINAL_285_TRUE_HOLDOUT_METRICS.csv` | Results | Holdout $C$-index and Brier score metrics | **APPROVED** | **NO** | `FINAL_285_TRUE_HOLDOUT_METRICS.csv` | Aggregate performance summary |
| `results/nested_cv_285_fold_metrics.csv` | Results | Fold-level CV metrics across 25 outer folds | **APPROVED** | **NO** | `nested_cv_285_fold_metrics.csv` | Aggregate CV performance summary |
| `results/FINAL_CLINICAL_BENCHMARK_METRICS.csv` | Results | Clinical vs Omics vs Combined benchmark summary | **APPROVED** | **NO** | `FINAL_CLINICAL_BENCHMARK_METRICS.csv` | Model comparison metrics |
| `results/FINAL_MOFA_FACTOR_SUMMARY.csv` | Results | Modality and factor variance decomposition | **APPROVED** | **NO** | `FINAL_MOFA_FACTOR_SUMMARY.csv` | $R^2$ variance explained per factor |
| `results/FINAL_F2_TOP_RNA_LOADINGS.csv` | Results | Top RNA loadings for Factor 2 | **APPROVED** | **NO** | `FINAL_F2_TOP_RNA_LOADINGS.csv` | Gene symbols and loading weights |
| `results/FINAL_FACTOR_CORRELATION_MATRIX_25RUNS.csv`| Results | Cross-fold factor alignment matrix | **APPROVED** | **NO** | `FINAL_FACTOR_CORRELATION_MATRIX_25RUNS.csv`| Correlation matrix |
| `results/FINAL_FACTOR_STABILITY_AUDIT.csv` | Results | Factor recovery and stability statistics | **APPROVED** | **NO** | `FINAL_FACTOR_STABILITY_AUDIT.csv` | Factor 2 stability ($|r| = 0.9784$) |
| `results/FINAL_PATHWAY_ENRICHMENT_AUDIT.csv` | Results | GO Biological Process and KEGG enrichment table | **APPROVED** | **NO** | `FINAL_PATHWAY_ENRICHMENT_AUDIT.csv` | Pathway names, $p$-values, $q$-values |
| `results/FINAL_BIOLOGICAL_CLAIM_CLASSIFICATION.csv` | Results | Curated biological claim classification | **APPROVED** | **NO** | `FINAL_BIOLOGICAL_CLAIM_CLASSIFICATION.csv`| Association vs causal claim audit |
| `results/nested_cv_285_projection_diagnostics.csv` | Results | Gram matrix condition numbers across folds | **APPROVED** | **NO** | `nested_cv_285_projection_diagnostics.csv`| Condition numbers |
| `code/01_preprocess_rna.R` | Code | RNA log-CPM and ComBat transformation | **APPROVED** | **NO** | `gate2b1_rna_pipeline.R` | Path-normalized to relative paths |
| `code/02_preprocess_methylation.R` | Code | Methylation probe filtering and ComBat | **APPROVED** | **NO** | `gate2b2_meth_pipeline.R` | Path-normalized to relative paths |
| `code/03_preprocess_cnv.R` | Code | Length-weighted segment-to-gene aggregation | **APPROVED** | **NO** | `gate2b3_cnv_pipeline.R` | Path-normalized; Affymetrix SNP 6.0 overlap |
| `code/04_preprocess_snv.R` | Code | Filtering to nine functional mutation classes | **APPROVED** | **NO** | `gate2b4_snv_pipeline.R` | Path-normalized |
| `code/05_generate_seed_fold_matrix.R` | Code | Deterministic 5x5 fold assignment generation | **APPROVED** | **NO** | `create_285_seed_fold_matrix.R` | Path-normalized |
| `code/06_run_nested_cv.R` | Code | Main 25-fold nested cross-validation engine | **APPROVED** | **NO** | `run_285_nested_cv.R` | Path-normalized; auto-detects Python |
| `code/07_train_and_evaluate_pipeline.R` | Code | Full Level 3 retraining and holdout evaluation pipeline | **APPROVED** | **NO** | `evaluate_285_production_and_holdout.R` | Requires upstream raw GDC inputs |
| `code/09_apply_frozen_model.R` | Code | Pure R Level 2 frozen model projection & prediction engine | **APPROVED** | **NO** | New Standalone Script | Zero Python dependency; zero raw data dependency |
| `code/08_compute_pooled_cindex.R` | Code | Cross-fold pooled Harrell's $C$-index calculation | **APPROVED** | **NO** | `compute_pooled_285_cindex.R` | Path-normalized |
| `code/build_master_manuscript.py` | Code | Markdown manuscript compiler | **APPROVED** | **NO** | `build_master_manuscript.py` | Path-normalized to relative paths |
| `code/verify_all_manuscript_facts.py` | Code | 32-metric and prohibited-claim red-team auditor | **APPROVED** | **NO** | `verify_all_manuscript_facts.py` | Path-normalized to relative paths |
| `figures/main/Figure_1` to `Figure_5` (PDF/PNG) | Figures | Publication main figures 1 through 5 | **APPROVED** | **NO** | `figures/main/` | 10 files (vector PDF and high-res PNG) |
| `figures/supplementary/Figure_S1` to `Figure_S8` | Figures | Publication supplementary figures S1 to S8 | **APPROVED** | **NO** | `figures/supplementary/` | 16 files (vector PDF and high-res PNG) |
| `docs/MANUSCRIPT_MASTER.md` | Documentation | Consolidated manuscript with certified audit | **APPROVED** | **NO** | `MANUSCRIPT_MASTER.md` | 163.5 KB complete master draft |
| `docs/MANUSCRIPT_EVIDENCE_MAP.md` | Documentation | Traceability matrix linking all claims to files | **APPROVED** | **NO** | `MANUSCRIPT_EVIDENCE_MAP.md` | Factual verification ledger |
| `docs/data_access.md` | Documentation | Instructions on GDC source data acquisition | **APPROVED** | **NO** | New Draft | Explains required GDC source inputs |
| `docs/reproducibility.md` | Documentation | Three-tiered reproducibility execution guide | **APPROVED** | **NO** | New Draft | Details Levels 1, 2, and 3 execution paths |
| `docs/manuscript_sections/` (01 to 07) | Documentation | Modular manuscript source files | **APPROVED** | **NO** | `manuscript_sections/` | 7 markdown section files |
