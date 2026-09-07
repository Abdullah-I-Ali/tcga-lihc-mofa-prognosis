# ==============================================================================
# NOTE: LEVEL 3 FULL 25-RUN NESTED CROSS-VALIDATION PIPELINE
# This script executes the complete 5-repeat 5-fold nested cross-validation.
# It requires preprocessed omics matrices from scripts 01-04 and survival vectors.
# To evaluate the pre-trained frozen model without retraining, use:
# code/09_apply_frozen_model.R
# ==============================================================================
###############################################################################
# run_285_nested_cv.R
# TCGA-LIHC MULTI-OMICS SURVIVAL PIPELINE
# CLEAN DEVELOPMENT COHORT (N = 285) 25-RUN NESTED CROSS-VALIDATION
#
# Cohort: N = 285 locked TCGA-LIHC primary tumor patients (103 Deaths, 182 Censored)
# Excludes the locked N=65 holdout completely (Intersection = 0)
#
# Architecture:
#   1. In-fold clamped ComBat (min_delta = 0.25)
#   2. In-fold feature selection on training data only
#   3. MOFA2 (15 factors) fit on outer training folds
#   4. Corrected Joint Ridge-Regularized Factor Projection (lambda = 1.0)
#   5. Models: ElasticNet-Cox (alpha = 0.5), Lasso-Cox (alpha = 1.0),
#      XGBoost-AFT, Random Survival Forest
#   6. Valid Breslow baseline hazard estimator & IPCW Brier/IBS calculations
#   7. Checkpointing per run into checkpoints_285_nested_cv/
###############################################################################

Sys.setenv(
  TMP = tempdir(),
  TEMP = tempdir(),
  RETICULATE_PYTHON = py_exe,
  OPENBLAS_NUM_THREADS = "4",
  OMP_NUM_THREADS = "4"
)

suppressPackageStartupMessages({
  library(data.table)
  library(caret)
  library(edgeR)
  library(matrixStats)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(reticulate)
  if (file.exists(py_exe)) use_python(py_exe, required = FALSE)
  library(MOFA2)
  library(survival)
  library(glmnet)
  library(xgboost)
  library(randomForestSRC)
})

cat("================================================================================\n")
cat("  CLEAN N = 285 DEVELOPMENT COHORT 25-RUN NESTED CROSS-VALIDATION\n")
cat("  Start Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================================\n\n")

# ------------------------------------------------------------------------------
# 1. Output Paths & Directory Setup (Strict Isolation from All Prior Runs)
# ------------------------------------------------------------------------------
chk_dir <- "checkpoints_285_nested_cv"
if (!dir.exists(chk_dir)) dir.create(chk_dir, recursive = TRUE)

manifest_file <- "nested_cv_285_execution_manifest.csv"
metrics_file  <- "nested_cv_285_fold_metrics.csv"
oof_file      <- "nested_cv_285_patient_oof_predictions.csv"
diag_file     <- "nested_cv_285_projection_diagnostics.csv"

if (!file.exists(manifest_file)) {
  fwrite(data.table(
    run_index = integer(), repeat_id = integer(), fold_id = integer(), run_id = character(),
    train_N = integer(), validation_N = integer(), death_train = integer(), death_validation = integer(),
    censored_train = integer(), censored_val = integer(),
    mofa_seed = integer(), inner_cv_seed = integer(), lasso_seed = integer(),
    xgboost_seed = integer(), rsf_seed = integer(), bootstrap_seed = integer(),
    status = character(), runtime_seconds = numeric(), mofa_time_seconds = numeric(),
    warnings_errors = character(), timestamp = character()
  ), manifest_file)
}

if (!file.exists(metrics_file)) {
  fwrite(data.table(
    run_index = integer(), run_id = character(), repeat_id = integer(), fold_id = integer(),
    model = character(), c_index = numeric(), 
    brier_1yr = numeric(), brier_3yr = numeric(), brier_5yr = numeric(), ibs = numeric()
  ), metrics_file)
}

if (!file.exists(oof_file)) {
  fwrite(data.table(
    patient_id = character(), run_index = integer(), repeat_id = integer(), fold_id = integer(), run_id = character(),
    os_time_days = numeric(), os_event = integer(),
    risk_lasso = numeric(), risk_enet = numeric(), risk_xgb = numeric(), risk_rf = numeric(),
    surv1yr_lasso = numeric(), surv3yr_lasso = numeric(), surv5yr_lasso = numeric(),
    surv1yr_enet = numeric(), surv3yr_enet = numeric(), surv5yr_enet = numeric(),
    surv1yr_xgb = numeric(), surv3yr_xgb = numeric(), surv5yr_xgb = numeric(),
    surv1yr_rf = numeric(), surv3yr_rf = numeric(), surv5yr_rf = numeric()
  ), oof_file)
}

if (!file.exists(diag_file)) {
  fwrite(data.table(
    run_index = integer(), run_id = character(), repeat_id = integer(), fold_id = integer(),
    omega_min_ev = numeric(), omega_max_ev = numeric(), omega_cond = numeric(),
    sd_ratio_mean = numeric(), sd_ratio_min = numeric(), sd_ratio_max = numeric(),
    z_val_min = numeric(), z_val_max = numeric(), nan_count = integer()
  ), diag_file)
}

# ------------------------------------------------------------------------------
# 2. Mathematical Helper Functions
# ------------------------------------------------------------------------------

fit_combat_train_clamped <- function(dat_tr, batch_tr, min_delta = 0.25, max_delta = 4.0) {
  batch_tr <- as.factor(batch_tr); batch_levels <- levels(batch_tr)
  grand_mean <- rowMeans(dat_tr); var_pooled <- rowVars(dat_tr)
  n_batch <- length(batch_levels)
  gamma_hat <- matrix(0, nrow(dat_tr), n_batch, dimnames = list(rownames(dat_tr), batch_levels))
  delta_hat <- matrix(1, nrow(dat_tr), n_batch, dimnames = list(rownames(dat_tr), batch_levels))
  for (b in batch_levels) {
    idx <- which(batch_tr == b); sub_dat <- dat_tr[, idx, drop = FALSE]
    gamma_hat[, b] <- rowMeans(sub_dat) - grand_mean
    if (length(idx) > 1) {
      d_hat <- sqrt(pmax(rowVars(sub_dat), 1e-8) / pmax(var_pooled, 1e-8))
      delta_hat[, b] <- pmin(pmax(d_hat, min_delta), max_delta)
    }
  }
  list(grand_mean = grand_mean, var_pooled = var_pooled, gamma_hat = gamma_hat,
       delta_hat = delta_hat, batch_levels = batch_levels)
}

apply_combat_frozen <- function(dat_new, batch_new, params) {
  dat_adj <- dat_new; batch_new <- as.character(batch_new)
  for (b in params$batch_levels) {
    idx <- which(batch_new == b)
    if (length(idx) > 0)
      dat_adj[, idx] <- (dat_new[, idx, drop = FALSE] - params$grand_mean -
                        params$gamma_hat[, b]) / params$delta_hat[, b] + params$grand_mean
  }
  dat_adj
}

calc_cindex <- function(surv_obj, risk_score) {
  as.numeric(concordance(surv_obj ~ I(-risk_score))$concordance)
}

# Valid Breslow Survival Probability Estimator for Cox Models
predict_cox_surv <- function(cv_fit, x_train, surv_train, x_val, eval_times = c(365, 1095, 1825)) {
  beta <- as.matrix(coef(cv_fit, s = "lambda.min"))
  eta_tr <- as.numeric(x_train %*% beta)
  eta_val <- as.numeric(x_val %*% beta)
  t_tr <- surv_train[, 1]; d_tr <- surv_train[, 2]
  event_times <- sort(unique(t_tr[d_tr == 1]))
  dH0 <- numeric(length(event_times))
  for (j in seq_along(event_times)) {
    tj <- event_times[j]
    dj <- sum(t_tr == tj & d_tr == 1)
    risk_sum <- sum(exp(eta_tr[t_tr >= tj]))
    dH0[j] <- dj / max(risk_sum, 1e-8)
  }
  H0 <- cumsum(dH0)
  h0_step <- stepfun(event_times, c(0, H0))
  surv_mat <- matrix(NA, nrow = nrow(x_val), ncol = length(eval_times))
  colnames(surv_mat) <- paste0("t_", eval_times)
  for (k in seq_along(eval_times)) {
    t_star <- eval_times[k]
    surv_mat[, k] <- exp(-h0_step(t_star) * exp(eta_val))
  }
  surv_mat <- pmin(pmax(surv_mat, 0), 1)
  for (k in 2:length(eval_times)) surv_mat[, k] <- pmin(surv_mat[, k], surv_mat[, k - 1])
  surv_mat
}

# Valid RSF Ensemble Survival Probability Curve
predict_rsf_surv <- function(rf_pred, eval_times = c(365, 1095, 1825)) {
  times_grid <- rf_pred$time.interest
  surv_curve <- rf_pred$survival
  surv_mat <- matrix(NA, nrow = nrow(surv_curve), ncol = length(eval_times))
  colnames(surv_mat) <- paste0("t_", eval_times)
  for (k in seq_along(eval_times)) {
    t_star <- eval_times[k]
    idx <- which(times_grid <= t_star)
    surv_mat[, k] <- if (length(idx) == 0) 1.0 else surv_curve[, max(idx)]
  }
  surv_mat <- pmin(pmax(surv_mat, 0), 1)
  for (k in 2:length(eval_times)) surv_mat[, k] <- pmin(surv_mat[, k], surv_mat[, k - 1])
  surv_mat
}

# Valid Parametric Normal AFT Survival Probability for XGBoost
predict_xgb_surv <- function(xgb_model, dval_aft, sigma = 1.20, eval_times = c(365, 1095, 1825)) {
  pred_mu <- as.numeric(stats::predict(xgb_model, dval_aft))
  surv_mat <- matrix(NA, nrow = length(pred_mu), ncol = length(eval_times))
  colnames(surv_mat) <- paste0("t_", eval_times)
  for (k in seq_along(eval_times)) {
    t_star <- eval_times[k]
    z_score <- (log(t_star) - pred_mu) / sigma
    surv_mat[, k] <- 1 - pnorm(z_score)
  }
  surv_mat <- pmin(pmax(surv_mat, 0), 1)
  for (k in 2:length(eval_times)) surv_mat[, k] <- pmin(surv_mat[, k], surv_mat[, k - 1])
  surv_mat
}

# Mathematically Valid IPCW Brier Score & IBS Calculation
calc_true_ipcw_brier <- function(train_surv, val_surv, pred_surv_mat, eval_times = c(365, 1095, 1825)) {
  c_fit <- survfit(Surv(train_surv[, 1], 1 - train_surv[, 2]) ~ 1)
  t_val <- val_surv[, 1]; d_val <- val_surv[, 2]
  bs_vec <- numeric(length(eval_times))
  names(bs_vec) <- paste0("BS_", eval_times)
  g_step <- stepfun(c_fit$time, c(1, c_fit$surv))
  for (k in seq_along(eval_times)) {
    t_star <- eval_times[k]
    g_ti <- pmax(g_step(t_val), 0.01)
    g_tstar <- max(g_step(t_star), 0.01)
    weights <- ifelse(t_val <= t_star & d_val == 1, 1 / g_ti,
               ifelse(t_val > t_star, 1 / g_tstar, 0))
    y_actual <- as.numeric(t_val > t_star)
    p_pred <- pred_surv_mat[, k]
    bs_vec[k] <- mean(weights * (y_actual - p_pred)^2)
  }
  list(brier_scores = bs_vec, ibs = mean(bs_vec))
}

# ------------------------------------------------------------------------------
# 3. Load Verified Clean Development Cohort (N = 285) & Raw Data
# ------------------------------------------------------------------------------
cat("Loading verified development cohort metadata (N = 285)...\n")
dev_ids_dt <- fread("./TRUE_DEV_285_IDS.csv")
locked_ids <- sort(dev_ids_dt$patient_id)
stopifnot(length(locked_ids) == 285)

# Verify zero holdout overlap
test_ids_quarantined <- fread("./FINAL_TEST_IDS.csv")$patient_id
stopifnot(length(intersect(locked_ids, test_ids_quarantined)) == 0)
cat(sprintf("Quarantined holdout verification: 0 overlap with %d holdout IDs.\n", length(test_ids_quarantined)))

event_vector <- fread("./TRUE_DEV_285_EVENT_VECTOR.csv")[match(locked_ids, patient_id)]
stopifnot(identical(event_vector$patient_id, locked_ids))
stopifnot(sum(event_vector$os_event == 1) == 103)
stopifnot(sum(event_vector$os_event == 0) == 182)
stopifnot(!any(is.na(event_vector$os_time_days)))

fold_matrix <- fread("./TRUE_DEV_285_SEED_FOLD_MATRIX.csv")

# Create reproducible outer multi-folds
set.seed(42)
multi_folds <- createMultiFolds(as.factor(event_vector$os_event), k = 5, times = 5)

# Ingest RNA
cat("Ingesting raw RNA-Seq data for N = 285...\n")
rna_raw <- readRDS("./rna_expression_raw.rds")
rna_mat_raw <- if (inherits(rna_raw, "SummarizedExperiment")) SummarizedExperiment::assay(rna_raw) else as.matrix(rna_raw)
storage.mode(rna_mat_raw) <- "double"
rna_manifest <- fread("./GATE_2B1_RNA_SAMPLE_MANIFEST.csv")
rna_aliqs <- rna_manifest[IN_LOCKED_PRIMARY == TRUE & PATIENT_ID %in% locked_ids, ALIQUOT_BARCODE]
names(rna_aliqs) <- rna_manifest[IN_LOCKED_PRIMARY == TRUE & PATIENT_ID %in% locked_ids, PATIENT_ID]
rna_mat_all <- rna_mat_raw[, rna_aliqs]
colnames(rna_mat_all) <- names(rna_aliqs)
rna_mat_all <- rna_mat_all[, locked_ids, drop = FALSE]
rm(rna_raw, rna_mat_raw); gc()

all_ens_clean <- sub("\\..*$", "", rownames(rna_mat_all))
gene_map_global <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(all_ens_clean), columns = "SYMBOL", keytype = "ENSEMBL")
gene_map_clean <- gene_map_global[!is.na(gene_map_global$SYMBOL) & gene_map_global$SYMBOL != "", ]
rm(gene_map_global); gc()

# Ingest Methylation
cat("Ingesting raw DNA Methylation 450k data for N = 285...\n")
meth_raw <- readRDS("./methylation_beta_raw.rds")
meth_mat_raw <- if (inherits(meth_raw, "SummarizedExperiment")) SummarizedExperiment::assay(meth_raw) else as.matrix(meth_raw)
storage.mode(meth_mat_raw) <- "double"
meth_manifest <- fread("./GATE_2B2_METH_SAMPLE_MANIFEST.csv")
meth_aliqs <- meth_manifest[IN_LOCKED_PRIMARY == TRUE & PATIENT_ID %in% locked_ids, ALIQUOT_BARCODE]
names(meth_aliqs) <- meth_manifest[IN_LOCKED_PRIMARY == TRUE & PATIENT_ID %in% locked_ids, PATIENT_ID]
meth_mat_all <- meth_mat_raw[, meth_aliqs]
colnames(meth_mat_all) <- names(meth_aliqs)
meth_mat_all <- meth_mat_all[, locked_ids, drop = FALSE]
rm(meth_raw, meth_mat_raw); gc()

# Apply static Chen 2013 + sex/control masking
chen_file <- "./Chen_2013_cross_reactive_probes.csv"
chen_probes <- if(file.exists(chen_file)) read.csv(chen_file, header = TRUE, stringsAsFactors = FALSE)[, 1] else character(0)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
chr_map <- anno$chr[match(rownames(meth_mat_all), rownames(anno))]
is_sex_chr     <- chr_map %in% c("chrX", "chrY")
is_non_cpg     <- grepl("^(rs|ch|ctl)", rownames(meth_mat_all), ignore.case = TRUE)
is_cross_react <- rownames(meth_mat_all) %in% chen_probes
static_mask_meth <- !is_cross_react & !is_sex_chr & !is_non_cpg
meth_mat_all <- meth_mat_all[static_mask_meth, , drop = FALSE]
rm(anno, chr_map, chen_probes); gc()

# Ingest CNV & SNV
cat("Ingesting CNV and SNV matrices for N = 285...\n")
cnv_full <- readRDS("./GATE_2B3_CNV_MATRIX.rds")$full_matrix[, locked_ids, drop = FALSE]
snv_full <- readRDS("./GATE_2B4_SNV_MATRIX.rds")$full_matrix[, locked_ids, drop = FALSE]

cat("All development data loaded successfully.\n\n")

# ------------------------------------------------------------------------------
# 4. Master 25-Run Execution Loop
# ------------------------------------------------------------------------------
total_runs <- 25
start_all_time <- proc.time()

for (r_idx in 1:total_runs) {
  run_start_time <- proc.time()
  row_info <- fold_matrix[run_index == r_idx]
  run_id_val <- row_info$run_id
  rep_id_val <- row_info$repeat_id
  fld_id_val <- row_info$fold_id
  
  chk_target <- file.path(chk_dir, sprintf("nested_cv_285_checkpoint_run_%d.rds", r_idx))
  
  cat(sprintf("\n================================================================================\n"))
  cat(sprintf(">>> STARTING RUN %d / %d: %s (Repeat %d, Fold %d)\n", 
              r_idx, total_runs, run_id_val, rep_id_val, fld_id_val))
  cat(sprintf("================================================================================\n"))
  
  # Resumability Checkpoint Verification
  if (file.exists(chk_target)) {
    cat(sprintf("Checkpoint %s already exists. Verifying integrity...\n", chk_target))
    chk_obj <- tryCatch(readRDS(chk_target), error = function(e) NULL)
    if (!is.null(chk_obj) && !is.null(chk_obj$metrics) && nrow(chk_obj$metrics) == 4) {
      cat(sprintf("Run %d verified intact. Skipping...\n", r_idx))
      next
    } else {
      cat(sprintf("Checkpoint corrupt or incomplete. Re-running run %d...\n", r_idx))
    }
  }
  
  # 1. Outer Train/Validation Split
  tr_idx  <- multi_folds[[r_idx]]
  val_idx <- setdiff(1:285, tr_idx)
  tr_ids  <- sort(locked_ids[tr_idx])
  val_ids <- sort(locked_ids[val_idx])
  
  stopifnot(length(intersect(tr_ids, val_ids)) == 0)
  stopifnot(length(union(tr_ids, val_ids)) == 285)
  
  ev_tr  <- event_vector[match(tr_ids, patient_id)]
  ev_val <- event_vector[match(val_ids, patient_id)]
  surv_tr  <- Surv(ev_tr$os_time_days, ev_tr$os_event)
  surv_val <- Surv(ev_val$os_time_days, ev_val$os_event)
  
  cat(sprintf("Split: Train N = %d (Deaths = %d, Censored = %d) | Val N = %d (Deaths = %d, Censored = %d)\n",
              length(tr_ids), sum(ev_tr$os_event == 1), sum(ev_tr$os_event == 0),
              length(val_ids), sum(ev_val$os_event == 1), sum(ev_val$os_event == 0)))
  
  # 2. RNA In-Fold Pipeline (Train Only Preprocessing)
  cat("Preprocessing RNA (in-fold)...\n")
  rna_raw_tr  <- rna_mat_all[, tr_ids, drop = FALSE]
  rna_raw_val <- rna_mat_all[, val_ids, drop = FALSE]
  
  non_zero_tr <- rowSums(rna_raw_tr) > 0
  rna_nz_tr <- rna_raw_tr[non_zero_tr, , drop = FALSE]
  cpm_raw_tr <- edgeR::cpm(rna_nz_tr)
  keep_expr_tr <- rowSums(cpm_raw_tr > 1.0) >= ceiling(0.20 * ncol(rna_nz_tr))
  rna_filt_tr <- rna_nz_tr[keep_expr_tr, , drop = FALSE]
  
  dge_tr <- DGEList(counts = rna_filt_tr)
  dge_tr <- calcNormFactors(dge_tr, method = "TMM")
  tr_norm_factors <- dge_tr$samples$norm.factors
  log2_cpm_tr <- edgeR::cpm(dge_tr, log = TRUE, prior.count = 2)
  
  clean_ens_tr <- sub("\\..*$", "", rownames(log2_cpm_tr))
  mapped_idx <- which(clean_ens_tr %in% gene_map_clean$ENSEMBL)
  log2_cpm_mapped <- log2_cpm_tr[mapped_idx, , drop = FALSE]
  clean_ens_mapped <- clean_ens_tr[mapped_idx]
  symbol_vec <- gene_map_clean$SYMBOL[match(clean_ens_mapped, gene_map_clean$ENSEMBL)]
  annot_df <- data.frame(ENSEMBL = clean_ens_mapped, ORIG_ENS = rownames(log2_cpm_mapped),
                         SYMBOL = symbol_vec, MEAN_EXPR = rowMeans(log2_cpm_mapped), stringsAsFactors = FALSE)
  annot_df <- annot_df[order(annot_df$SYMBOL, -annot_df$MEAN_EXPR), ]
  annot_unique <- annot_df[!duplicated(annot_df$SYMBOL), ]
  log2_cpm_sym_tr <- log2_cpm_mapped[annot_unique$ORIG_ENS, , drop = FALSE]
  rownames(log2_cpm_sym_tr) <- annot_unique$SYMBOL
  
  batch_tr_rna <- substr(rna_aliqs[colnames(log2_cpm_sym_tr)], 22, 25)
  combat_params_rna <- fit_combat_train_clamped(log2_cpm_sym_tr, batch_tr_rna)
  rna_tr_combat <- apply_combat_frozen(log2_cpm_sym_tr, batch_tr_rna, combat_params_rna)
  
  active_genes_tr <- which(rowMeans(rna_tr_combat) >= 1.0)
  top_rna_idx <- active_genes_tr[order(rowMads(rna_tr_combat)[active_genes_tr],
                                       decreasing = TRUE)[seq_len(min(902, length(active_genes_tr)))]]
  rna_sel_genes <- sort(rownames(rna_tr_combat)[top_rna_idx])
  rna_train_mat <- rna_tr_combat[rna_sel_genes, , drop = FALSE]
  
  # Project RNA Val with Frozen Parameters
  raw_val_mapped <- rna_raw_val[annot_unique$ORIG_ENS, , drop = FALSE]
  rownames(raw_val_mapped) <- annot_unique$SYMBOL
  lib_sizes_val <- colSums(rna_raw_val)
  eff_lib_val <- lib_sizes_val * median(tr_norm_factors)
  cpm_val_raw <- t(t((raw_val_mapped + 2) / (eff_lib_val + 2)) * 1e6)
  log2_cpm_val <- log2(cpm_val_raw)
  batch_val_rna <- substr(rna_aliqs[colnames(log2_cpm_val)], 22, 25)
  rna_val_mat <- apply_combat_frozen(log2_cpm_val, batch_val_rna, combat_params_rna)[rna_sel_genes, , drop = FALSE]
  
  # 3. Methylation In-Fold Pipeline (Optimized Method B)
  cat("Preprocessing Methylation (in-fold Method B)...\n")
  meth_raw_tr  <- meth_mat_all[, tr_ids, drop = FALSE]
  meth_raw_val <- meth_mat_all[, val_ids, drop = FALSE]
  
  tr_miss_frac <- rowMeans(is.na(meth_raw_tr))
  keep_probes_m <- rownames(meth_raw_tr)[tr_miss_frac <= 0.05]
  meth_tr_sub  <- meth_raw_tr[keep_probes_m, , drop = FALSE]
  meth_val_sub <- meth_raw_val[keep_probes_m, , drop = FALSE]
  rm(meth_raw_tr, meth_raw_val); gc()
  
  tr_meth_meds <- rowMedians(meth_tr_sub, na.rm = TRUE)
  na_idx_tr <- which(is.na(meth_tr_sub), arr.ind = TRUE)
  if (nrow(na_idx_tr) > 0) meth_tr_sub[na_idx_tr] <- tr_meth_meds[na_idx_tr[, 1]]
  na_idx_val <- which(is.na(meth_val_sub), arr.ind = TRUE)
  if (nrow(na_idx_val) > 0) meth_val_sub[na_idx_val] <- tr_meth_meds[na_idx_val[, 1]]
  
  clamp_beta <- function(b) pmin(pmax(b, 1e-4), 1 - 1e-4)
  beta2M     <- function(b) { bc <- clamp_beta(b); log2(bc / (1 - bc)) }
  m_tr  <- beta2M(meth_tr_sub)
  m_val <- beta2M(meth_val_sub)
  rm(meth_tr_sub, meth_val_sub, na_idx_tr, na_idx_val); gc()
  
  probe_vars_tr <- rowVars(m_tr)
  n_meth_sel <- min(5000, length(probe_vars_tr))
  meth_sel_probes <- sort(rownames(m_tr)[order(probe_vars_tr, decreasing = TRUE)[seq_len(n_meth_sel)]])
  m_tr_sel  <- m_tr[meth_sel_probes, , drop = FALSE]
  m_val_sel <- m_val[meth_sel_probes, , drop = FALSE]
  rm(m_tr, m_val); gc()
  
  batch_tr_meth <- substr(meth_aliqs[colnames(m_tr_sel)], 22, 25)
  batch_val_meth <- substr(meth_aliqs[colnames(m_val_sel)], 22, 25)
  combat_params_meth <- fit_combat_train_clamped(m_tr_sel, batch_tr_meth)
  meth_train_mat <- apply_combat_frozen(m_tr_sel, batch_tr_meth, combat_params_meth)
  meth_val_mat   <- apply_combat_frozen(m_val_sel, batch_val_meth, combat_params_meth)
  rm(m_tr_sel, m_val_sel); gc()
  
  # 4. CNV & SNV In-Fold Feature Selection
  cnv_tr  <- cnv_full[, tr_ids, drop = FALSE]
  cnv_val <- cnv_full[, val_ids, drop = FALSE]
  cnv_sel_genes <- sort(rownames(cnv_tr)[order(rowVars(cnv_tr), decreasing = TRUE)[seq_len(min(5000, nrow(cnv_tr)))]])
  cnv_train_mat <- cnv_tr[cnv_sel_genes, , drop = FALSE]
  cnv_val_mat   <- cnv_val[cnv_sel_genes, , drop = FALSE]
  
  snv_tr  <- snv_full[, tr_ids, drop = FALSE]
  snv_val <- snv_full[, val_ids, drop = FALSE]
  snv_sel_genes <- sort(rownames(snv_tr)[order(rowSums(snv_tr), decreasing = TRUE)[seq_len(min(500, nrow(snv_tr)))]])
  snv_train_mat <- snv_tr[snv_sel_genes, , drop = FALSE]
  snv_val_mat   <- snv_val[snv_sel_genes, , drop = FALSE]
  
  # 5. Fit MOFA2 on Outer Train (15 Factors)
  cat("Fitting MOFA2 on training data (15 factors)...\n")
  mofa_views <- list(RNA = rna_train_mat, Methylation = meth_train_mat, CNV = cnv_train_mat, SNV = snv_train_mat)
  mofa_obj <- create_mofa(mofa_views)
  
  data_opts <- get_default_data_options(mofa_obj)
  data_opts$scale_views <- FALSE; data_opts$scale_groups <- FALSE
  model_opts <- get_default_model_options(mofa_obj)
  model_opts$num_factors <- 15
  model_opts$likelihoods['RNA']         <- 'gaussian'
  model_opts$likelihoods['Methylation'] <- 'gaussian'
  model_opts$likelihoods['CNV']         <- 'gaussian'
  model_opts$likelihoods['SNV']         <- 'bernoulli'
  
  train_opts <- get_default_training_options(mofa_obj)
  train_opts$seed             <- row_info$mofa_seed
  train_opts$maxiter          <- 1000L
  train_opts$convergence_mode <- 'fast'
  train_opts$verbose          <- FALSE
  
  mofa_obj <- prepare_mofa(mofa_obj, data_options = data_opts, model_options = model_opts, training_options = train_opts)
  
  t0_m <- proc.time()
  mofa_fit <- run_mofa(mofa_obj, use_basilisk = FALSE)
  mofa_dur <- (proc.time() - t0_m)[["elapsed"]]
  
  Z_train <- get_factors(mofa_fit)[[1]]
  W_list  <- get_weights(mofa_fit)
  
  # 6. Corrected Joint Ridge-Regularized Factor Projection (lambda = 1.0)
  cat("Computing joint ridge-regularized factor projection (lambda = 1.0)...\n")
  U_train  <- matrix(0, nrow = ncol(rna_train_mat), ncol = 15)
  U_val    <- matrix(0, nrow = ncol(rna_val_mat),   ncol = 15)
  Omega_tr <- matrix(0, 15, 15)
  
  for (vname in c("RNA", "Methylation", "CNV", "SNV")) {
    X_tr_v  <- switch(vname, RNA = rna_train_mat, Methylation = meth_train_mat, CNV = cnv_train_mat, SNV = snv_train_mat)
    X_val_v <- switch(vname, RNA = rna_val_mat,   Methylation = meth_val_mat,   CNV = cnv_val_mat,   SNV = snv_val_mat)
    W_v     <- W_list[[vname]]
    
    clean_w_names <- sub("_(RNA|Methylation|CNV|SNV)$", "", rownames(W_v))
    common_f <- intersect(clean_w_names, rownames(X_tr_v))
    
    w_idx <- match(common_f, clean_w_names)
    x_idx <- match(common_f, rownames(X_tr_v))
    
    W_sub     <- W_v[w_idx, , drop = FALSE]
    X_tr_sub  <- X_tr_v[x_idx, , drop = FALSE]
    X_val_sub <- X_val_v[x_idx, , drop = FALSE]
    
    mu_tr_v    <- rowMeans(X_tr_sub)
    X_tr_cent  <- X_tr_sub  - mu_tr_v
    X_val_cent <- X_val_sub - mu_tr_v
    
    U_train  <- U_train  + t(X_tr_cent)  %*% W_sub
    U_val    <- U_val    + t(X_val_cent) %*% W_sub
    Omega_tr <- Omega_tr + crossprod(W_sub)
  }
  
  lambda_ridge <- 1.0
  inv_omega <- solve(Omega_tr + lambda_ridge * diag(15))
  Z_val <- U_val %*% inv_omega
  
  colnames(Z_train) <- paste0("Factor", 1:ncol(Z_train))
  colnames(Z_val)   <- paste0("Factor", 1:ncol(Z_val))
  rownames(Z_train) <- tr_ids
  rownames(Z_val)   <- val_ids
  
  # Diagnostics
  ev_omega <- eigen(Omega_tr)$values
  sd_ratios <- apply(Z_val, 2, sd) / apply(Z_train, 2, sd)
  diag_row <- data.table(
    run_index = r_idx, run_id = run_id_val, repeat_id = rep_id_val, fold_id = fld_id_val,
    omega_min_ev = min(ev_omega), omega_max_ev = max(ev_omega), omega_cond = max(ev_omega)/min(ev_omega),
    sd_ratio_mean = mean(sd_ratios), sd_ratio_min = min(sd_ratios), sd_ratio_max = max(sd_ratios),
    z_val_min = min(Z_val), z_val_max = max(Z_val), nan_count = sum(!is.finite(Z_val))
  )
  fwrite(diag_row, diag_file, append = TRUE)
  
  cat(sprintf("Diagnostics: Cond = %.1f | SD ratio = %.3f [%.3f, %.3f] | Range = [%.2f, %.2f]\n",
              diag_row$omega_cond, diag_row$sd_ratio_mean, diag_row$sd_ratio_min, diag_row$sd_ratio_max,
              diag_row$z_val_min, diag_row$z_val_max))
  
  stopifnot(diag_row$nan_count == 0)
  stopifnot(diag_row$sd_ratio_mean < 10.0)
  
  # 7. Model 1: LASSO-Cox
  cat("Fitting LASSO-Cox...\n")
  set.seed(row_info$lasso_seed)
  cv_lasso <- cv.glmnet(x = Z_train, y = surv_tr, family = "cox", alpha = 1, nfolds = 5)
  risk_lasso <- as.numeric(stats::predict(cv_lasso, newx = Z_val, s = "lambda.min", type = "link"))
  c_lasso <- calc_cindex(surv_val, risk_lasso)
  surv_p_lasso <- predict_cox_surv(cv_lasso, Z_train, surv_tr, Z_val)
  brier_lasso <- calc_true_ipcw_brier(surv_tr, surv_val, surv_p_lasso)
  
  # 8. Model 2: ElasticNet-Cox (Primary Model)
  cat("Fitting ElasticNet-Cox (alpha = 0.5)...\n")
  set.seed(row_info$lasso_seed)
  cv_enet <- cv.glmnet(x = Z_train, y = surv_tr, family = "cox", alpha = 0.5, nfolds = 5)
  risk_enet <- as.numeric(stats::predict(cv_enet, newx = Z_val, s = "lambda.min", type = "link"))
  c_enet <- calc_cindex(surv_val, risk_enet)
  surv_p_enet <- predict_cox_surv(cv_enet, Z_train, surv_tr, Z_val)
  brier_enet <- calc_true_ipcw_brier(surv_tr, surv_val, surv_p_enet)
  
  # 9. Model 3: XGBoost-AFT
  cat("Fitting XGBoost-AFT...\n")
  set.seed(row_info$xgboost_seed)
  dtrain_aft <- xgb.DMatrix(data = Z_train)
  setinfo(dtrain_aft, "label_lower_bound", ev_tr$os_time_days)
  setinfo(dtrain_aft, "label_upper_bound", ifelse(ev_tr$os_event == 1, ev_tr$os_time_days, Inf))
  dval_aft <- xgb.DMatrix(data = Z_val)
  xgb_params_aft <- list(
    objective = "survival:aft", eval_metric = "aft-nloglik",
    aft_loss_distribution = "normal", aft_loss_distribution_scale = 1.20,
    max_depth = 3, eta = 0.05, subsample = 0.8
  )
  cv_xgb <- xgb.cv(params = xgb_params_aft, data = dtrain_aft, nrounds = 100, nfold = 5,
                   early_stopping_rounds = 10, verbose = 0)
  best_nrounds <- cv_xgb$best_iteration
  if (is.null(best_nrounds) || best_nrounds < 1) best_nrounds <- 20
  xgb_model <- xgb.train(params = xgb_params_aft, data = dtrain_aft, nrounds = best_nrounds, verbose = 0)
  risk_xgb <- -as.numeric(stats::predict(xgb_model, dval_aft))
  c_xgb <- calc_cindex(surv_val, risk_xgb)
  surv_p_xgb <- predict_xgb_surv(xgb_model, dval_aft, sigma = 1.20)
  brier_xgb <- calc_true_ipcw_brier(surv_tr, surv_val, surv_p_xgb)
  
  # 10. Model 4: Random Survival Forest
  cat("Fitting Random Survival Forest...\n")
  set.seed(row_info$rsf_seed)
  df_train <- as.data.frame(Z_train)
  df_train$os_time  <- ev_tr$os_time_days
  df_train$os_event <- ev_tr$os_event
  df_val <- as.data.frame(Z_val)
  rf_fit <- rfsrc(Surv(os_time, os_event) ~ ., data = df_train, ntree = 200, splitrule = "logrank", seed = row_info$rsf_seed)
  rf_pred <- stats::predict(rf_fit, newdata = df_val)
  risk_rf <- rf_pred$predicted
  c_rf <- calc_cindex(surv_val, risk_rf)
  surv_p_rf <- predict_rsf_surv(rf_pred)
  brier_rf <- calc_true_ipcw_brier(surv_tr, surv_val, surv_p_rf)
  
  # Assemble Out-of-Fold Records
  oof_dt <- data.table(
    patient_id = val_ids, run_index = r_idx, repeat_id = rep_id_val, fold_id = fld_id_val, run_id = run_id_val,
    os_time_days = ev_val$os_time_days, os_event = ev_val$os_event,
    risk_lasso = risk_lasso, risk_enet = risk_enet, risk_xgb = risk_xgb, risk_rf = risk_rf,
    surv1yr_lasso = surv_p_lasso[, 1], surv3yr_lasso = surv_p_lasso[, 2], surv5yr_lasso = surv_p_lasso[, 3],
    surv1yr_enet  = surv_p_enet[, 1],  surv3yr_enet  = surv_p_enet[, 2],  surv5yr_enet  = surv_p_enet[, 3],
    surv1yr_xgb   = surv_p_xgb[, 1],   surv3yr_xgb   = surv_p_xgb[, 2],   surv5yr_xgb   = surv_p_xgb[, 3],
    surv1yr_rf    = surv_p_rf[, 1],    surv3yr_rf    = surv_p_rf[, 2],    surv5yr_rf    = surv_p_rf[, 3]
  )
  fwrite(oof_dt, oof_file, append = TRUE)
  
  # Assemble Metrics Records
  metrics_dt <- data.table(
    run_index = r_idx, run_id = run_id_val, repeat_id = rep_id_val, fold_id = fld_id_val,
    model = c("LASSO-Cox", "ElasticNet", "XGBoost-AFT", "RandomSurvivalForest"),
    c_index = c(c_lasso, c_enet, c_xgb, c_rf),
    brier_1yr = c(brier_lasso$brier_scores[1], brier_enet$brier_scores[1], brier_xgb$brier_scores[1], brier_rf$brier_scores[1]),
    brier_3yr = c(brier_lasso$brier_scores[2], brier_enet$brier_scores[2], brier_xgb$brier_scores[2], brier_rf$brier_scores[2]),
    brier_5yr = c(brier_lasso$brier_scores[3], brier_enet$brier_scores[3], brier_xgb$brier_scores[3], brier_rf$brier_scores[3]),
    ibs       = c(brier_lasso$ibs,             brier_enet$ibs,             brier_xgb$ibs,             brier_rf$ibs)
  )
  fwrite(metrics_dt, metrics_file, append = TRUE)
  
  run_dur <- (proc.time() - run_start_time)[["elapsed"]]
  
  # Assemble Manifest Record
  man_row <- data.table(
    run_index = r_idx, repeat_id = rep_id_val, fold_id = fld_id_val, run_id = run_id_val,
    train_N = length(tr_ids), validation_N = length(val_ids),
    death_train = sum(ev_tr$os_event == 1), death_validation = sum(ev_val$os_event == 1),
    censored_train = sum(ev_tr$os_event == 0), censored_val = sum(ev_val$os_event == 0),
    mofa_seed = row_info$mofa_seed, inner_cv_seed = row_info$inner_cv_seed, lasso_seed = row_info$lasso_seed,
    xgboost_seed = row_info$xgboost_seed, rsf_seed = row_info$rsf_seed, bootstrap_seed = row_info$bootstrap_seed,
    status = "SUCCESS", runtime_seconds = round(run_dur, 2), mofa_time_seconds = round(mofa_dur, 2),
    warnings_errors = "NONE", timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  fwrite(man_row, manifest_file, append = TRUE)
  
  checkpoint_data <- list(
    manifest_row = man_row,
    oof_predictions = oof_dt,
    metrics = metrics_dt,
    diagnostics = diag_row,
    selected_features = list(RNA = rna_sel_genes, Methylation = meth_sel_probes, CNV = cnv_sel_genes, SNV = snv_sel_genes),
    Z_train = Z_train,
    Z_val = Z_val,
    Omega_tr = Omega_tr,
    inv_omega = inv_omega
  )
  saveRDS(checkpoint_data, chk_target)
  
  cat(sprintf(">>> COMPLETED RUN %d in %.1f s (MOFA: %.1f s) | C-indices: LASSO=%.4f, ENet=%.4f, XGB=%.4f, RSF=%.4f | IBS: LASSO=%.4f, ENet=%.4f\n",
              r_idx, run_dur, mofa_dur, c_lasso, c_enet, c_xgb, c_rf, brier_lasso$ibs, brier_enet$ibs))
  
  rm(rna_train_mat, rna_val_mat, meth_train_mat, meth_val_mat, cnv_train_mat, cnv_val_mat, snv_train_mat, snv_val_mat,
     mofa_obj, mofa_fit, Z_train, Z_val, W_list, U_train, U_val, Omega_tr, inv_omega,
     cv_lasso, cv_enet, xgb_model, rf_fit, oof_dt, metrics_dt, checkpoint_data)
  gc()
}

total_elapsed <- (proc.time() - start_all_time)[["elapsed"]]
cat("\n================================================================================\n")
cat(sprintf("ALL 25 RUNS ON N = 285 COMPLETE in %.1f minutes!\n", total_elapsed / 60))
cat("================================================================================\n")

all_metrics <- fread(metrics_file)
all_oof     <- fread(oof_file)
all_diag    <- fread(diag_file)
all_man     <- fread(manifest_file)

consolidated_results <- list(
  manifest = all_man,
  fold_metrics = all_metrics,
  oof_predictions = all_oof,
  diagnostics = all_diag,
  total_runtime_minutes = total_elapsed / 60,
  completion_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
saveRDS(consolidated_results, "nested_cv_285_results.rds")
cat("Consolidated results saved to nested_cv_285_results.rds\n")
