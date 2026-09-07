# Data Access and Raw TCGA-LIHC Source Acquisition

This repository provides the code, sanitized model parameters, feature definitions, and evaluation workflows for the study:
**"A Reproducibility-Focused Multi-Omics Latent-Factor Framework for Prognostic Modeling in Hepatocellular Carcinoma"**

---

## 1. Public Availability of Source Data

All primary molecular and clinical data analyzed in this study derive from **The Cancer Genome Atlas Liver Hepatocellular Carcinoma (TCGA-LIHC)** project and are hosted publicly by the National Cancer Institute (NCI) Genomic Data Commons (GDC):
- **GDC Project Portal:** [https://portal.gdc.cancer.gov/projects/TCGA-LIHC](https://portal.gdc.cancer.gov/projects/TCGA-LIHC)
- **Primary Open-Access Accession:** `TCGA-LIHC`

The patient-level raw sequencing reads, methylation arrays, copy number segmentations, somatic mutation calls, and clinical records are subject to standard GDC open-access data use terms and are not redistributed directly within this repository.

---

## 2. Required Input Datasets for Full Upstream Preprocessing

To execute the raw feature preprocessing pipelines (`code/01_preprocess_rna.R` through `code/04_preprocess_snv.R`), the following upstream files must be acquired from GDC:

| Modality | Data Type / Format | GDC Workflow / Data Category | Expected Local Input File |
| :--- | :--- | :--- | :--- |
| **Gene Expression (RNA)** | HTSeq - Counts or STAR - Counts (Raw Read Counts) | Transcriptome Profiling | `rna_expression_raw.rds` |
| **DNA Methylation** | Infinium HumanMethylation450 Beta Values / IDATs | DNA Methylation | `methylation_beta_raw.rds` |
| **Copy Number Variation** | Affymetrix Genome-Wide Human SNP Array 6.0 Segment Mean | Copy Number Variation | `cnv_segment_raw.rds` |
| **Somatic Mutations (SNV)** | Somatic MAF (MuTect2 Variant Calls) | Sequencing Reads / Simple Nucleotide Variation | `snv_mutation_raw.rds` |
| **Clinical Metadata** | TCGA-LIHC Clinical Supplement & TCGA-CDR Survival Data | Clinical | `clinical_data.tsv` |

---

## 3. Current Documentation Status & Planned Manifest Release

> [!IMPORTANT]
> **A direct GDC Data Transfer Tool manifest (`gdc_manifest.txt`) is pending final release upon preprint publication.**
> 
> Currently, the repository provides:
> 1. Exact cohort patient identifiers in `manifests/TRUE_DEV_285_IDS.csv` ($N = 285$) and `manifests/FINAL_TEST_IDS.csv` ($N = 65$).
> 2. Reference gene coordinates (`data/references/gene_coords_hg38.rds`).
> 3. Cross-reactive methylation probe filters (`data/references/Chen_2013_cross_reactive_probes.csv`).
> 4. Deterministic cross-validation fold allocations (`data/references/TRUE_DEV_285_SEED_FOLD_MATRIX.csv`).
> 5. The sanitized, frozen multi-omics production model (`model/FINAL_285_FROZEN_MODEL_BUNDLE_SANITIZED.rds`).
>
> Full end-to-end raw-data download automation will be finalized with a locked GDC UUID manifest alongside the preprint release. In the interim, downstream model evaluation, out-of-sample regularized projection math, and statistical benchmark verification can be executed directly using the provided model bundle and summary result tables.
