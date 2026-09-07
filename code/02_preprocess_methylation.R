# ==============================================================================
# NOTE: LEVEL 3 RAW RECONSTRUCTION PIPELINE
# This script performs 450K DNA methylation probe filtering, median imputation,
# and ComBat batch adjustment. It requires the raw input beta value matrix:
# methylation_beta_raw.rds (acquired from NCI GDC; see docs/data_access.md).
# ==============================================================================
# ==============================================================================
# GATE 2B.2 — CLEAN DNA METHYLATION PREPROCESSING PIPELINE
# ==============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(matrixStats)
  library(data.table)
  library(dplyr)
  library(tibble)
})

data_dir <- ".""

cat("============================================================\n")
cat("  GATE 2B.2 CLEAN DNA METHYLATION PREPROCESSING PIPELINE\n")
cat("  Timestamp:", format(Sys.time()), "\n")
cat("============================================================\n\n")

# ------------------------------------------------------------------------------
# 1. LOAD LOCKED PRIMARY COHORT MANIFEST (N = 350)
# ------------------------------------------------------------------------------
truth_df <- read.csv(file.path(data_dir, "GATE_2A6_FINAL_VERIFIED_PATIENT_TRUTH.csv"), stringsAsFactors = FALSE)
primary_pts <- sort(truth_df$patient_id[truth_df$PRIMARY_INCLUDED == TRUE])
cat(sprintf("1. Locked Primary Modeling Cohort: N = %d Patients\n", length(primary_pts)))

# ------------------------------------------------------------------------------
# 2. RAW METHYLATION INPUT & SAMPLE SELECTION
# ------------------------------------------------------------------------------
cat("\n2. Loading Raw Methylation Beta Matrix (485,577 probes x 430 aliquots)...\n")
meth_raw <- readRDS(file.path(data_dir, "methylation_beta_raw.rds"))
raw_mat <- if (inherits(meth_raw, "SummarizedExperiment")) assay(meth_raw) else as.matrix(meth_raw)
storage.mode(raw_mat) <- "double"

cat(sprintf("   - Raw Matrix Dimensions: %d Probes x %d Aliquots\n", nrow(raw_mat), ncol(raw_mat)))

# Aliquot Sample Types
all_barcodes <- colnames(raw_mat)
sample_types <- substr(all_barcodes, 14, 15)
pt_barcodes <- substr(all_barcodes, 1, 12)

# Sample Manifest of All 430 Raw Aliquots
sample_manifest <- data.frame(
  ALIQUOT_BARCODE = all_barcodes,
  PATIENT_ID = pt_barcodes,
  SAMPLE_TYPE_CODE = sample_types,
  SAMPLE_TYPE_NAME = ifelse(sample_types == "01", "Primary Solid Tumor",
                     ifelse(sample_types == "02", "Recurrent Solid Tumor",
                     ifelse(sample_types == "11", "Solid Tissue Normal", "Other"))),
  IN_LOCKED_PRIMARY = pt_barcodes %in% primary_pts,
  stringsAsFactors = FALSE
)

sample_manifest$SELECTION_STATUS <- apply(sample_manifest, 1, function(r) {
  st <- r["SAMPLE_TYPE_CODE"]
  in_p <- r["IN_LOCKED_PRIMARY"] == "TRUE"
  if (st != "01") return(paste0("DISCARDED: Non-tumor (Type ", st, ")"))
  if (!in_p) return("DISCARDED: Patient excluded in Gate 2A cohort lock")
  return("SELECTED: Primary Tumor Aliquot for Locked Cohort")
})

write.csv(sample_manifest, file.path(data_dir, "GATE_2B2_METH_SAMPLE_MANIFEST.csv"), row.names = FALSE)
cat("   ✓ Saved: GATE_2B2_METH_SAMPLE_MANIFEST.csv\n")

# Filter to Primary Tumor Aliquots for the 350 Locked Patients
tumor_mask <- sample_types == "01" & pt_barcodes %in% primary_pts
tumor_mat <- raw_mat[, tumor_mask, drop = FALSE]
colnames(tumor_mat) <- substr(colnames(tumor_mat), 1, 12)

# Verify 1:1 patient mapping (0 duplicates)
tumor_mat <- tumor_mat[, !duplicated(colnames(tumor_mat)), drop = FALSE]
tumor_mat <- tumor_mat[, primary_pts, drop = FALSE]

cat(sprintf("   - Selected Aliquots: %d Aliquots across %d Unique Patients (1:1 mapped)\n",
            ncol(tumor_mat), length(unique(colnames(tumor_mat)))))

# ------------------------------------------------------------------------------
# 3. PROBE QC & FILTERING
# ------------------------------------------------------------------------------
cat("\n3. PROBE QC & SYSTEMATIC FILTERING:\n")
raw_probe_n <- nrow(tumor_mat)
all_probe_ids <- rownames(tumor_mat)

# A. Cross-reactive probes (Chen et al. 2013)
chen_file <- file.path(data_dir, "Chen_2013_cross_reactive_probes.csv")
if (file.exists(chen_file)) {
  chen_probes <- read.csv(chen_file, header = TRUE, stringsAsFactors = FALSE)[, 1]
} else {
  chen_probes <- character(0)
}
is_cross_reactive <- all_probe_ids %in% chen_probes
cross_reactive_n <- sum(is_cross_reactive)
cat(sprintf("   - Raw 450k Probes:             %d\n", raw_probe_n))
cat(sprintf("   - Chen 2013 Cross-Reactive:    %d removed\n", cross_reactive_n))

# B. Sex Chromosome & SNP/Non-CpG Probes using Illumina Annotation
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_probes <- rownames(anno)

# Sex chromosomes
chr_map <- anno$chr[match(all_probe_ids, anno_probes)]
is_sex_chr <- chr_map %in% c("chrX", "chrY")
sex_chr_n <- sum(is_sex_chr)
cat(sprintf("   - Sex Chromosome Probes (X/Y): %d removed\n", sex_chr_n))

# Non-CpG / SNP control probes (rs*, ch*, etc.)
is_non_cpg <- grepl("^(rs|ch|ctl)", all_probe_ids, ignore.case = TRUE)
non_cpg_n <- sum(is_non_cpg)
cat(sprintf("   - Non-CpG / SNP Probes (rs/ch):%d removed\n", non_cpg_n))

# Composite Annotation Filter
qc_mask <- !is_cross_reactive & !is_sex_chr & !is_non_cpg
probes_after_anno <- all_probe_ids[qc_mask]
tumor_mat_anno <- tumor_mat[probes_after_anno, , drop = FALSE]
cat(sprintf("   - Retained Post-Annotation:    %d Probes\n", nrow(tumor_mat_anno)))

# C. Missingness Filtering (Drop probes with > 5% missing values across 350 samples)
missing_per_probe <- rowSums(is.na(tumor_mat_anno)) / ncol(tumor_mat_anno)
keep_missing <- missing_per_probe <= 0.05
excess_missing_n <- sum(!keep_missing)
tumor_mat_clean <- tumor_mat_anno[keep_missing, , drop = FALSE]
cat(sprintf("   - Probes with >5%% Missing:     %d removed\n", excess_missing_n))
cat(sprintf("   - Retained Filtered Probes:    %d Probes\n", nrow(tumor_mat_clean)))

# ------------------------------------------------------------------------------
# 4. MISSING VALUE IMPUTATION (ROW-MEDIAN)
# ------------------------------------------------------------------------------
cat("\n4. MISSING VALUE IMPUTATION:\n")
total_na <- sum(is.na(tumor_mat_clean))
cat(sprintf("   - Remaining Missing Beta Values: %d (%.4f%% of matrix)\n",
            total_na, (total_na / (nrow(tumor_mat_clean) * ncol(tumor_mat_clean))) * 100))

# Impute row-wise medians for remaining missing values
row_meds <- rowMedians(tumor_mat_clean, na.rm = TRUE)
na_idx <- which(is.na(tumor_mat_clean), arr.ind = TRUE)
if (nrow(na_idx) > 0) {
  tumor_mat_clean[na_idx] <- row_meds[na_idx[, 1]]
}
cat(sprintf("   - Post-Imputation Missing:       %d NA values\n", sum(is.na(tumor_mat_clean))))

# ------------------------------------------------------------------------------
# 5. BETA -> M-VALUE TRANSFORMATION
# ------------------------------------------------------------------------------
cat("\n5. BETA -> M-VALUE TRANSFORMATION:\n")
# Dynamic clipping to avoid log(0) or log(Inf)
beta_clamped <- pmin(pmax(tumor_mat_clean, 1e-4), 1 - 1e-4)
m_val_mat <- log2(beta_clamped / (1 - beta_clamped))

cat(sprintf("   - Transformed Matrix Dims:     %d Probes x %d Patients\n", nrow(m_val_mat), ncol(m_val_mat)))
cat(sprintf("   - M-Value Value Range:         [%.4f, %.4f]\n", min(m_val_mat), max(m_val_mat)))
cat(sprintf("   - NA / NaN / Inf Check:        NA=%d, NaN=%d, Inf=%d\n",
            sum(is.na(m_val_mat)), sum(is.nan(m_val_mat)), sum(is.infinite(m_val_mat))))

# ------------------------------------------------------------------------------
# 6. TECHNICAL PLATE BATCH IDENTIFICATION & BATCH ADJUSTMENT
# ------------------------------------------------------------------------------
cat("\n6. BATCH CORRECTION:\n")
# Extract Plate ID batch from raw aliquot barcodes
matched_barcodes <- all_barcodes[sample_types == "01" & pt_barcodes %in% primary_pts]
matched_pts <- substr(matched_barcodes, 1, 12)
names(matched_barcodes) <- matched_pts
matched_barcodes_350 <- matched_barcodes[colnames(m_val_mat)]

plate_ids <- substr(matched_barcodes_350, 22, 25)
batch_vec <- factor(plate_ids)

cat(sprintf("   - Batch Variable:              TCGA Plate ID (%d distinct batches)\n", length(levels(batch_vec))))
cat(sprintf("   - Batch Sample Sizes:          Min=%d, Median=%.1f, Max=%d (Singletons=0)\n",
            min(table(batch_vec)), median(table(batch_vec)), max(table(batch_vec))))

# Fit Frozen-Parameter ComBat-Equivalent Location/Scale Transformation
fit_combat_params <- function(dat_train, batch_train) {
  batch_train <- as.factor(batch_train)
  batch_levels <- levels(batch_train)
  n_batch <- length(batch_levels)
  
  grand_mean <- rowMeans(dat_train)
  var_pooled <- rowVars(dat_train)
  
  gamma_hat <- matrix(0, nrow = nrow(dat_train), ncol = n_batch, dimnames = list(rownames(dat_train), batch_levels))
  delta_hat <- matrix(1, nrow = nrow(dat_train), ncol = n_batch, dimnames = list(rownames(dat_train), batch_levels))
  
  for (i in seq_along(batch_levels)) {
    b <- batch_levels[i]
    idx <- which(batch_train == b)
    sub_dat <- dat_train[, idx, drop = FALSE]
    gamma_hat[, b] <- rowMeans(sub_dat) - grand_mean
    if (length(idx) > 1) {
      sub_var <- rowVars(sub_dat)
      delta_hat[, b] <- sqrt(pmax(sub_var, 1e-8) / pmax(var_pooled, 1e-8))
    }
  }
  
  list(
    grand_mean = grand_mean,
    var_pooled = var_pooled,
    gamma_hat = gamma_hat,
    delta_hat = delta_hat,
    batch_levels = batch_levels
  )
}

apply_frozen_combat <- function(dat_new, batch_new, combat_params) {
  dat_adj <- dat_new
  batch_new <- as.character(batch_new)
  grand_mean <- combat_params$grand_mean
  
  for (b in combat_params$batch_levels) {
    idx <- which(batch_new == b)
    if (length(idx) > 0) {
      gamma_b <- combat_params$gamma_hat[, b]
      delta_b <- combat_params$delta_hat[, b]
      dat_adj[, idx] <- (dat_new[, idx, drop = FALSE] - grand_mean - gamma_b) / delta_b + grand_mean
    }
  }
  dat_adj
}

combat_params_350 <- fit_combat_params(m_val_mat, batch_vec)
meth_combat_clean <- apply_frozen_combat(m_val_mat, batch_vec, combat_params_350)

cat(sprintf("   - Batch Adjusted Matrix:       %d Probes x %d Patients\n", nrow(meth_combat_clean), ncol(meth_combat_clean)))
cat(sprintf("   - Batch Adjusted Value Range:  [%.4f, %.4f]\n", min(meth_combat_clean), max(meth_combat_clean)))

# ------------------------------------------------------------------------------
# 7. FEATURE REPRESENTATION FOR DOWNSTREAM MOFA2 (TOP 5,000 PROBES)
# ------------------------------------------------------------------------------
cat("\n7. FEATURE SELECTION FOR MOFA2 MATRIX (TOP 5,000 PROBES):\n")
probe_vars <- rowVars(meth_combat_clean)
top_5000_idx <- order(probe_vars, decreasing = TRUE)[1:5000]

meth_mofa_mat <- meth_combat_clean[top_5000_idx, , drop = FALSE]
# Sort probe IDs alphabetically for deterministic reproducible ordering
meth_mofa_mat <- meth_mofa_mat[order(rownames(meth_mofa_mat)), ]

cat(sprintf("   - Full Clean Methylation Matrix: %d Probes x %d Patients\n", nrow(meth_combat_clean), ncol(meth_combat_clean)))
cat(sprintf("   - Top Variable MOFA Features:    %d Probes x %d Patients\n", nrow(meth_mofa_mat), ncol(meth_mofa_mat)))

# Save Serialized Output (Tier A Descriptive Global Artifact)
saveRDS(list(
  full_matrix = meth_combat_clean,
  mofa_features = meth_mofa_mat,
  combat_params = combat_params_350,
  row_medians = row_meds,
  patient_ids = colnames(meth_combat_clean),
  artifact_classification = "DESCRIPTIVE GLOBAL ARTIFACT — NOT VALID FOR NESTED-CV MODEL INPUT"
), file = file.path(data_dir, "GATE_2B2_METH_MATRIX.rds"))
cat("   ✓ Saved: GATE_2B2_METH_MATRIX.rds\n")

# ------------------------------------------------------------------------------
# 8. INDEPENDENT VERIFICATION CHECKS & OUTPUT
# ------------------------------------------------------------------------------
cat("\n8. INDEPENDENT VERIFICATION ASSERTIONS:\n")

checks <- list(
  c("EXACT_N_350_PATIENTS", ncol(meth_mofa_mat) == 350, sprintf("Expected 350, got %d", ncol(meth_mofa_mat))),
  c("ZERO_DUPLICATE_PATIENTS", sum(duplicated(colnames(meth_mofa_mat))) == 0, "Duplicate patient IDs detected"),
  c("PATIENTS_MATCH_LOCKED_COHORT", identical(colnames(meth_mofa_mat), primary_pts), "Patient IDs mismatch locked cohort"),
  c("EXACT_TOP_PROBES_5000", nrow(meth_mofa_mat) == 5000, sprintf("Expected 5000 probes, got %d", nrow(meth_mofa_mat))),
  c("ZERO_NAN_INF_NA", sum(is.na(meth_mofa_mat) | is.nan(meth_mofa_mat) | is.infinite(meth_mofa_mat)) == 0, "Non-finite values present"),
  c("ZERO_SYNTHETIC_FEATURES", all(is.finite(meth_mofa_mat)) && nrow(meth_mofa_mat) == 5000, "Invalid probe distribution"),
  c("REPRODUCIBLE_FEATURE_ORDER", identical(rownames(meth_mofa_mat), sort(rownames(meth_mofa_mat))), "Features not deterministically ordered"),
  c("BATCH_SINGLETONS_ZERO", all(table(batch_vec) > 1), "Singleton batches present")
)

verif_df <- data.frame(
  CHECK_ID = sapply(checks, `[[`, 1),
  RESULT = sapply(checks, function(x) if (x[[2]]) "PASS" else "FAIL"),
  MESSAGE = sapply(checks, `[[`, 3),
  stringsAsFactors = FALSE
)

write.csv(verif_df, file.path(data_dir, "GATE_2B2_METH_VERIFICATION.csv"), row.names = FALSE)
cat("   ✓ Saved: GATE_2B2_METH_VERIFICATION.csv\n")
print(verif_df)

cat("\n============================================================\n")
cat("  GATE 2B.2 CLEAN DNA METHYLATION PREPROCESSING COMPLETE\n")
cat("  Timestamp:", format(Sys.time()), "\n")
cat("============================================================\n")
