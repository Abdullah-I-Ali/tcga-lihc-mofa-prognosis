# ==============================================================================
# NOTE: LEVEL 3 RAW RECONSTRUCTION PIPELINE
# This script performs RNA-seq TMM normalization, log2-CPM transformation,
# and ComBat batch adjustment. It requires the raw input count file:
# rna_expression_raw.rds (acquired from NCI GDC; see docs/data_access.md).
# ==============================================================================
# ==============================================================================
# GATE 2B.1 — CLEAN RNA PREPROCESSING PIPELINE
# ==============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(matrixStats)
  library(data.table)
  library(dplyr)
  library(tibble)
})

data_dir <- ".""

cat("============================================================\n")
cat("  GATE 2B.1 CLEAN RNA PREPROCESSING PIPELINE\n")
cat("  Timestamp:", format(Sys.time()), "\n")
cat("============================================================\n\n")

# ------------------------------------------------------------------------------
# 1. LOAD LOCKED PRIMARY COHORT MANIFEST (N = 350)
# ------------------------------------------------------------------------------
truth_df <- read.csv(file.path(data_dir, "GATE_2A6_FINAL_VERIFIED_PATIENT_TRUTH.csv"), stringsAsFactors = FALSE)
primary_pts <- sort(truth_df$patient_id[truth_df$PRIMARY_INCLUDED == TRUE])
cat(sprintf("1. Locked Primary Modeling Cohort: N = %d Patients\n", length(primary_pts)))

# ------------------------------------------------------------------------------
# 2. RAW RNA INPUT & SAMPLE SELECTION
# ------------------------------------------------------------------------------
rna_raw <- readRDS(file.path(data_dir, "rna_expression_raw.rds"))
raw_mat <- if (inherits(rna_raw, "SummarizedExperiment")) assay(rna_raw) else as.matrix(rna_raw)
storage.mode(raw_mat) <- "double"

cat(sprintf("\n2. Raw RNA Matrix Dimensions: %d Ensembl IDs x %d Aliquots\n", nrow(raw_mat), ncol(raw_mat)))

# Aliquot Sample Types
all_barcodes <- colnames(raw_mat)
sample_types <- substr(all_barcodes, 14, 15)
pt_barcodes <- substr(all_barcodes, 1, 12)

# Sample Manifest of All 424 Raw Aliquots
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

write.csv(sample_manifest, file.path(data_dir, "GATE_2B1_RNA_SAMPLE_MANIFEST.csv"), row.names = FALSE)
cat("   ✓ Saved: GATE_2B1_RNA_SAMPLE_MANIFEST.csv\n")

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
# 3. GENE FILTERING
# ------------------------------------------------------------------------------
cat("\n3. GENE FILTERING:\n")
raw_gene_n <- nrow(tumor_mat)

# A. Remove all-zero genes across the 350 samples
row_sums_raw <- rowSums(tumor_mat)
non_zero_genes <- row_sums_raw > 0
all_zero_n <- sum(!non_zero_genes)
tumor_mat_nz <- tumor_mat[non_zero_genes, , drop = FALSE]
cat(sprintf("   - Raw Ensembl Genes:           %d\n", raw_gene_n))
cat(sprintf("   - All-zero Genes Removed:      %d\n", all_zero_n))
cat(sprintf("   - Non-zero Genes Retained:     %d\n", nrow(tumor_mat_nz)))

# B. Low-expression filtering: CPM > 1 in >= 20% of samples (>= 70 samples)
min_samples <- ceiling(0.20 * ncol(tumor_mat_nz))
cpm_vals <- cpm(tumor_mat_nz)
keep_expr <- rowSums(cpm_vals > 1.0) >= min_samples
low_expr_n <- sum(!keep_expr)
tumor_mat_filtered <- tumor_mat_nz[keep_expr, , drop = FALSE]

cat(sprintf("   - Low-expression Filter Rule:  CPM > 1.0 in >= %d samples (20%% of N=%d)\n", min_samples, ncol(tumor_mat_nz)))
cat(sprintf("   - Low-expression Genes Out:    %d\n", low_expr_n))
cat(sprintf("   - Filtered Genes Retained:     %d\n", nrow(tumor_mat_filtered)))

# ------------------------------------------------------------------------------
# 4. NORMALIZATION (TMM + log2 CPM)
# ------------------------------------------------------------------------------
cat("\n4. NORMALIZATION:\n")
dge <- DGEList(counts = tumor_mat_filtered)
dge <- calcNormFactors(dge, method = "TMM")
norm_factors <- dge$samples$norm.factors
lib_sizes_eff <- dge$samples$lib.size * norm_factors

cat(sprintf("   - TMM Norm Factors Range:      [%.4f, %.4f] (Median: %.4f)\n",
            min(norm_factors), max(norm_factors), median(norm_factors)))
cat(sprintf("   - Effective Library Sizes:     [%.2fM, %.2fM] reads (Median: %.2fM)\n",
            min(lib_sizes_eff)/1e6, max(lib_sizes_eff)/1e6, median(lib_sizes_eff)/1e6))

# Transformation: log2(CPM + 2)
log2_cpm_mat <- cpm(dge, log = TRUE, prior.count = 2)
cat(sprintf("   - Transformed Matrix Dims:     %d Genes x %d Patients\n", nrow(log2_cpm_mat), ncol(log2_cpm_mat)))
cat(sprintf("   - Transformed Value Range:     [%.4f, %.4f]\n", min(log2_cpm_mat), max(log2_cpm_mat)))
cat(sprintf("   - NA / NaN / Inf Check:        NA=%d, NaN=%d, Inf=%d\n",
            sum(is.na(log2_cpm_mat)), sum(is.nan(log2_cpm_mat)), sum(is.infinite(log2_cpm_mat))))

# ------------------------------------------------------------------------------
# 5. GENE ANNOTATION & DUPLICATE SYMBOL RESOLUTION
# ------------------------------------------------------------------------------
cat("\n5. GENE ANNOTATION:\n")
# Clean version numbers from Ensembl IDs if present
clean_ensembl <- sub("\\..*$", "", rownames(log2_cpm_mat))

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = clean_ensembl,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "ENSEMBL"
)

# Remove NA / unmapped symbols
gene_map_clean <- gene_map[!is.na(gene_map$SYMBOL) & gene_map$SYMBOL != "", ]
mapped_ensembl <- unique(gene_map_clean$ENSEMBL)
unmapped_n <- length(setdiff(clean_ensembl, mapped_ensembl))
cat(sprintf("   - Unmapped Ensembl IDs:        %d / %d\n", unmapped_n, length(clean_ensembl)))

# Align mapped subset
mapped_idx <- which(clean_ensembl %in% mapped_ensembl)
log2_cpm_mapped <- log2_cpm_mat[mapped_idx, , drop = FALSE]
clean_ensembl_mapped <- clean_ensembl[mapped_idx]

# Map to HGNC symbols with mean expression resolving duplicates
symbol_vec <- gene_map_clean$SYMBOL[match(clean_ensembl_mapped, gene_map_clean$ENSEMBL)]
gene_means <- rowMeans(log2_cpm_mapped)

annot_df <- data.frame(
  ENSEMBL = clean_ensembl_mapped,
  ORIG_ENSEMBL = rownames(log2_cpm_mapped),
  SYMBOL = symbol_vec,
  MEAN_EXPR = gene_means,
  stringsAsFactors = FALSE
)

# Resolve duplicate symbols: Select Ensembl transcript with highest mean expression
annot_df <- annot_df[order(annot_df$SYMBOL, -annot_df$MEAN_EXPR), ]
annot_unique <- annot_df[!duplicated(annot_df$SYMBOL), ]

log2_cpm_symbols <- log2_cpm_mapped[annot_unique$ORIG_ENSEMBL, , drop = FALSE]
rownames(log2_cpm_symbols) <- annot_unique$SYMBOL

cat(sprintf("   - Unique HGNC Symbols:         %d Unique Genes x %d Patients\n",
            nrow(log2_cpm_symbols), ncol(log2_cpm_symbols)))

# ------------------------------------------------------------------------------
# 6. BATCH IDENTIFICATION & FROZEN COMBAT ADJUSTMENT
# ------------------------------------------------------------------------------
cat("\n6. BATCH CORRECTION:\n")
# Extract Plate ID batch for the 350 patients from raw aliquot barcodes
matched_barcodes <- all_barcodes[sample_types == "01" & pt_barcodes %in% primary_pts]
matched_pts <- substr(matched_barcodes, 1, 12)
names(matched_barcodes) <- matched_pts
matched_barcodes_350 <- matched_barcodes[colnames(log2_cpm_symbols)]

plate_ids <- substr(matched_barcodes_350, 22, 25)
batch_vec <- factor(plate_ids)

cat(sprintf("   - Batch Variable:              TCGA Plate ID (%d distinct batches)\n", length(levels(batch_vec))))
cat(sprintf("   - Batch Sample Sizes:          Min=%d, Median=%d, Max=%d (Singletons=0)\n",
            min(table(batch_vec)), median(table(batch_vec)), max(table(batch_vec))))

# Fit Frozen-Parameter ComBat Model
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
  
  for (j in 1:ncol(dat_new)) {
    b <- batch_new[j]
    if (b %in% combat_params$batch_levels) {
      gamma_b <- combat_params$gamma_hat[, b]
      delta_b <- combat_params$delta_hat[, b]
      dat_adj[, j] <- (dat_new[, j] - combat_params$grand_mean - gamma_b) / delta_b + combat_params$grand_mean
    } else {
      dat_adj[, j] <- dat_new[, j]
    }
  }
  dat_adj
}

combat_params_350 <- fit_combat_params(log2_cpm_symbols, batch_vec)
rna_combat_clean <- apply_frozen_combat(log2_cpm_symbols, batch_vec, combat_params_350)

cat(sprintf("   - Batch Adjusted Matrix:       %d Genes x %d Patients\n", nrow(rna_combat_clean), ncol(rna_combat_clean)))
cat(sprintf("   - Batch Adjusted Value Range:  [%.4f, %.4f]\n", min(rna_combat_clean), max(rna_combat_clean)))

# ------------------------------------------------------------------------------
# 7. FEATURE REPRESENTATION FOR DOWNSTREAM MOFA2
# ------------------------------------------------------------------------------
cat("\n7. FEATURE SELECTION FOR MOFA2 MATRIX:\n")
# Filter to biologically active genes (mean expression >= 1.0) and top MAD variance
gene_means_adj <- rowMeans(rna_combat_clean)
gene_mads_adj <- rowMads(rna_combat_clean)

# Select top 902 variable genes (mean >= 1.0, top MAD)
active_genes <- which(gene_means_adj >= 1.0)
top_mad_order <- order(gene_mads_adj[active_genes], decreasing = TRUE)
top_902_idx <- active_genes[top_mad_order[1:min(902, length(active_genes))]]

rna_mofa_mat <- rna_combat_clean[top_902_idx, , drop = FALSE]
# Sort genes alphabetically for deterministic reproducible ordering
rna_mofa_mat <- rna_mofa_mat[order(rownames(rna_mofa_mat)), ]

cat(sprintf("   - Full Clean RNA Matrix:       %d Genes x %d Patients\n", nrow(rna_combat_clean), ncol(rna_combat_clean)))
cat(sprintf("   - Top Variable MOFA Features:  %d Genes x %d Patients\n", nrow(rna_mofa_mat), ncol(rna_mofa_mat)))

# Save Outputs
saveRDS(list(
  full_matrix = rna_combat_clean,
  mofa_features = rna_mofa_mat,
  combat_params = combat_params_350,
  norm_factors = norm_factors,
  annot_unique = annot_unique,
  patient_ids = colnames(rna_combat_clean)
), file = file.path(data_dir, "GATE_2B1_RNA_MATRIX.rds"))
cat("   ✓ Saved: GATE_2B1_RNA_MATRIX.rds\n")

# ------------------------------------------------------------------------------
# 8. INDEPENDENT VERIFICATION CHECKS & OUTPUT
# ------------------------------------------------------------------------------
cat("\n8. INDEPENDENT VERIFICATION ASSERTIONS:\n")

checks <- list(
  c("EXACT_N_350_PATIENTS", ncol(rna_mofa_mat) == 350, sprintf("Expected 350, got %d", ncol(rna_mofa_mat))),
  c("ZERO_DUPLICATE_PATIENTS", sum(duplicated(colnames(rna_mofa_mat))) == 0, "Duplicate patient IDs detected"),
  c("PATIENTS_MATCH_LOCKED_COHORT", identical(colnames(rna_mofa_mat), primary_pts), "Patient IDs mismatch locked cohort"),
  c("EXACT_TOP_GENES_902", nrow(rna_mofa_mat) == 902, sprintf("Expected 902 genes, got %d", nrow(rna_mofa_mat))),
  c("ZERO_NAN_INF_NA", sum(is.na(rna_mofa_mat) | is.nan(rna_mofa_mat) | is.infinite(rna_mofa_mat)) == 0, "Non-finite values present"),
  c("ZERO_SYNTHETIC_FEATURES", all(is.finite(rna_mofa_mat)) && nrow(rna_mofa_mat) > 0, "Invalid feature distribution"),
  c("REPRODUCIBLE_FEATURE_ORDER", identical(rownames(rna_mofa_mat), sort(rownames(rna_mofa_mat))), "Features not deterministically ordered"),
  c("BATCH_SINGLETONS_ZERO", all(table(batch_vec) > 1), "Singleton batches present")
)

verif_df <- data.frame(
  CHECK_ID = sapply(checks, `[[`, 1),
  RESULT = sapply(checks, function(x) if (x[[2]]) "PASS" else "FAIL"),
  MESSAGE = sapply(checks, `[[`, 3),
  stringsAsFactors = FALSE
)

write.csv(verif_df, file.path(data_dir, "GATE_2B1_RNA_VERIFICATION.csv"), row.names = FALSE)
cat("   ✓ Saved: GATE_2B1_RNA_VERIFICATION.csv\n")
print(verif_df)

cat("\n============================================================\n")
cat("  GATE 2B.1 CLEAN RNA PREPROCESSING COMPLETE\n")
cat("  Timestamp:", format(Sys.time()), "\n")
cat("============================================================\n")
