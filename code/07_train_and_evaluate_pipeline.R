# ==============================================================================
# NOTE: LEVEL 3 FULL RECONSTRUCTION & TRAINING PIPELINE
# This script represents the complete training, factorization, and holdout
# evaluation pipeline from raw molecular inputs. It requires raw TCGA-LIHC inputs
# (rna_expression_raw.rds, methylation_beta_raw.rds, cnv_segment_raw.rds,
# snv_mutation_raw.rds) which must be acquired from GDC (see docs/data_access.md).
# To apply the pre-trained frozen model without raw data, use:
# code/09_apply_frozen_model.R
# ==============================================================================
###############################################################################
# evaluate_285_production_and_holdout.R
# TCGA-LIHC CLEAN PRODUCTION FIT (N = 285) & TRUE INDEPENDENT HOLDOUT EVAL (N = 65)
#
# Steps:
#   1. Verify completion of 25-run nested CV on N = 285
#   2. Fit final frozen production pipeline on ALL 285 development patients:
#      - In-cohort clamped ComBat & feature selection
#      - MOFA2 (15 factors) on N = 285
#      - Compute and freeze projection operator: (Omega_tr + 1.0 * I)^(-1)
#      - ElasticNet-Cox (alpha = 0.5) fit on Z_train with 5-fold CV to pick lambda.min
#      - Breslow baseline cumulative hazard H0(t) from N = 285
#      - Save bundle: FINAL_285_FROZEN_MODEL_BUNDLE.rds
#   3. Unlock and evaluate quarantined N = 65 holdout ONCE:
#      - Load FINAL_TEST_IDS.csv (N = 65, 19 deaths, 46 censored)
#      - Apply frozen N = 285 preprocessing & projection parameters
#      - Generate predictions: FINAL_285_TRUE_HOLDOUT_PREDICTIONS.csv
#      - Calculate metrics, 10,000 bootstrap CI, valid IPCW Brier/IBS, KM HR & log-rank
#      - Save metrics: FINAL_285_TRUE_HOLDOUT_METRICS.csv
###############################################################################

# Detect Python executable from environment or virtualenv
py_exe <- Sys.getenv("RETICULATE_PYTHON")
if (py_exe == "" || !file.exists(py_exe)) py_exe <- Sys.which("python")
if (py_exe == "" || !file.exists(py_exe)) py_exe <- Sys.which("python3")

Sys.setenv(
  TMP = tempdir(),
  TEMP = tempdir(),
  RETICULATE_PYTHON = py_exe,
  OPENBLAS_NUM_THREADS = "4",
  OMP_NUM_THREADS = "4"
)

suppressPackageStartupMessages({
  library(data.table)
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
})

cat("================================================================================\n")
cat("  TCGA-LIHC PRODUCTION REFIT (N = 285) & TRUE HOLDOUT EVALUATION (N = 65)\n")
cat("  Start Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================================\n\n")

# Helper functions
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

# 1. Load Dev IDs and Holdout IDs
dev_ids_dt <- fread("./TRUE_DEV_285_IDS.csv")
dev_ids <- sort(dev_ids_dt$patient_id)
stopifnot(length(dev_ids) == 285)

test_ids_dt <- fread("./FINAL_TEST_IDS.csv")
test_ids <- sort(test_ids_dt$patient_id)
stopifnot(length(test_ids) == 65)
stopifnot(length(intersect(dev_ids, test_ids)) == 0)

dev_events <- fread("./TRUE_DEV_285_EVENT_VECTOR.csv")[match(dev_ids, patient_id)]
test_events <- test_ids_dt[match(test_ids, patient_id)]

surv_train <- Surv(dev_events$os_time_days, dev_events$os_event)
surv_test  <- Surv(test_events$os_time, test_events$os_event)

# 2. Ingest Multi-Omics Data for all 350 patients (Dev + Test)
cat("Loading multi-omics inputs...\n")
rna_raw <- readRDS("./rna_expression_raw.rds")
rna_mat_raw <- if (inherits(rna_raw, "SummarizedExperiment")) SummarizedExperiment::assay(rna_raw) else as.matrix(rna_raw)
storage.mode(rna_mat_raw) <- "double"
rna_manifest <- fread("./GATE_2B1_RNA_SAMPLE_MANIFEST.csv")
all_350_ids <- sort(c(dev_ids, test_ids))
rna_aliqs <- rna_manifest[IN_LOCKED_PRIMARY == TRUE & PATIENT_ID %in% all_350_ids, ALIQUOT_BARCODE]
names(rna_aliqs) <- rna_manifest[IN_LOCKED_PRIMARY == TRUE & PATIENT_ID %in% all_350_ids, PATIENT_ID]
rna_mat_all <- rna_mat_raw[, rna_aliqs]
colnames(rna_mat_all) <- names(rna_aliqs)
rm(rna_raw, rna_mat_raw); gc()

all_ens_clean <- sub("\\..*$", "", rownames(rna_mat_all))
gene_map_global <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(all_ens_clean), columns = "SYMBOL", keytype = "ENSEMBL")
gene_map_clean <- gene_map_global[!is.na(gene_map_global$SYMBOL) & gene_map_global$SYMBOL != "", ]
rm(gene_map_global); gc()

meth_raw <- readRDS("./methylation_beta_raw.rds")
meth_mat_raw <- if (inherits(meth_raw, "SummarizedExperiment")) SummarizedExperiment::assay(meth_raw) else as.matrix(meth_raw)
storage.mode(meth_mat_raw) <- "double"
meth_manifest <- fread("./GATE_2B2_METH_SAMPLE_MANIFEST.csv")
meth_aliqs <- meth_manifest[IN_LOCKED_PRIMARY == TRUE & PATIENT_ID %in% all_350_ids, ALIQUOT_BARCODE]
names(meth_aliqs) <- meth_manifest[IN_LOCKED_PRIMARY == TRUE & PATIENT_ID %in% all_350_ids, PATIENT_ID]
meth_mat_all <- meth_mat_raw[, meth_aliqs]
colnames(meth_mat_all) <- names(meth_aliqs)
rm(meth_raw, meth_mat_raw); gc()

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

cnv_full <- readRDS("./GATE_2B3_CNV_MATRIX.rds")$full_matrix
snv_full <- readRDS("./GATE_2B4_SNV_MATRIX.rds")$full_matrix

# 3. Fit Production Pipeline on N = 285 ONLY
cat("\n>>> FITTING PRODUCTION PIPELINE ON N = 285 DEVELOPMENT COHORT <<<\n")

# RNA Preprocessing on Dev
rna_raw_tr <- rna_mat_all[, dev_ids, drop = FALSE]
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

# Methylation Preprocessing on Dev
meth_raw_tr <- meth_mat_all[, dev_ids, drop = FALSE]
tr_miss_frac <- rowMeans(is.na(meth_raw_tr))
keep_probes_m <- rownames(meth_raw_tr)[tr_miss_frac <= 0.05]
meth_tr_sub  <- meth_raw_tr[keep_probes_m, , drop = FALSE]
rm(meth_raw_tr); gc()

tr_meth_meds <- rowMedians(meth_tr_sub, na.rm = TRUE)
na_idx_tr <- which(is.na(meth_tr_sub), arr.ind = TRUE)
if (nrow(na_idx_tr) > 0) meth_tr_sub[na_idx_tr] <- tr_meth_meds[na_idx_tr[, 1]]

clamp_beta <- function(b) pmin(pmax(b, 1e-4), 1 - 1e-4)
beta2M     <- function(b) { bc <- clamp_beta(b); log2(bc / (1 - bc)) }
m_tr <- beta2M(meth_tr_sub)
rm(meth_tr_sub, na_idx_tr); gc()

probe_vars_tr <- rowVars(m_tr)
n_meth_sel <- min(5000, length(probe_vars_tr))
meth_sel_probes <- sort(rownames(m_tr)[order(probe_vars_tr, decreasing = TRUE)[seq_len(n_meth_sel)]])
m_tr_sel <- m_tr[meth_sel_probes, , drop = FALSE]
rm(m_tr); gc()

batch_tr_meth <- substr(meth_aliqs[colnames(m_tr_sel)], 22, 25)
combat_params_meth <- fit_combat_train_clamped(m_tr_sel, batch_tr_meth)
meth_train_mat <- apply_combat_frozen(m_tr_sel, batch_tr_meth, combat_params_meth)
rm(m_tr_sel); gc()

# CNV & SNV Feature Selection on Dev
cnv_tr <- cnv_full[, dev_ids, drop = FALSE]
cnv_sel_genes <- sort(rownames(cnv_tr)[order(rowVars(cnv_tr), decreasing = TRUE)[seq_len(min(5000, nrow(cnv_tr)))]])
cnv_train_mat <- cnv_tr[cnv_sel_genes, , drop = FALSE]

snv_tr <- snv_full[, dev_ids, drop = FALSE]
snv_sel_genes <- sort(rownames(snv_tr)[order(rowSums(snv_tr), decreasing = TRUE)[seq_len(min(500, nrow(snv_tr)))]])
snv_train_mat <- snv_tr[snv_sel_genes, , drop = FALSE]

# Fit Production MOFA2 (15 factors) on Dev
cat("Fitting production MOFA2 on N = 285...\n")
mofa_views_tr <- list(RNA = rna_train_mat, Methylation = meth_train_mat, CNV = cnv_train_mat, SNV = snv_train_mat)
mofa_obj_tr <- create_mofa(mofa_views_tr)

data_opts <- get_default_data_options(mofa_obj_tr)
data_opts$scale_views <- FALSE; data_opts$scale_groups <- FALSE
model_opts <- get_default_model_options(mofa_obj_tr)
model_opts$num_factors <- 15
model_opts$likelihoods['RNA']         <- 'gaussian'
model_opts$likelihoods['Methylation'] <- 'gaussian'
model_opts$likelihoods['CNV']         <- 'gaussian'
model_opts$likelihoods['SNV']         <- 'bernoulli'

train_opts <- get_default_training_options(mofa_obj_tr)
train_opts$seed             <- 42001
train_opts$maxiter          <- 1000L
train_opts$convergence_mode <- 'fast'
train_opts$verbose          <- FALSE

mofa_obj_tr <- prepare_mofa(mofa_obj_tr, data_options = data_opts, model_options = model_opts, training_options = train_opts)
mofa_fit_tr <- run_mofa(mofa_obj_tr, use_basilisk = FALSE)

Z_train <- get_factors(mofa_fit_tr)[[1]]
W_list  <- get_weights(mofa_fit_tr)
colnames(Z_train) <- paste0("Factor", 1:15)
rownames(Z_train) <- dev_ids

# Calculate & Freeze Projection Operator
Omega_tr <- matrix(0, 15, 15)
mu_tr_list <- list()
W_sub_list <- list()

for (vname in c("RNA", "Methylation", "CNV", "SNV")) {
  X_tr_v <- switch(vname, RNA = rna_train_mat, Methylation = meth_train_mat, CNV = cnv_train_mat, SNV = snv_train_mat)
  W_v    <- W_list[[vname]]
  clean_w <- sub("_(RNA|Methylation|CNV|SNV)$", "", rownames(W_v))
  common_f <- intersect(clean_w, rownames(X_tr_v))
  
  w_idx <- match(common_f, clean_w)
  x_idx <- match(common_f, rownames(X_tr_v))
  
  W_sub <- W_v[w_idx, , drop = FALSE]
  X_tr_sub <- X_tr_v[x_idx, , drop = FALSE]
  
  mu_tr_v <- rowMeans(X_tr_sub)
  mu_tr_list[[vname]] <- mu_tr_v
  W_sub_list[[vname]] <- W_sub
  
  Omega_tr <- Omega_tr + crossprod(W_sub)
}

lambda_ridge <- 1.0
inv_omega <- solve(Omega_tr + lambda_ridge * diag(15))

# Fit Production ElasticNet-Cox (alpha = 0.5) on Dev
cat("Fitting production ElasticNet-Cox (alpha = 0.5) on N = 285...\n")
set.seed(42002)
cv_enet_prod <- cv.glmnet(x = Z_train, y = surv_train, family = "cox", alpha = 0.5, nfolds = 5)
lambda_star <- cv_enet_prod$lambda.min
beta_enet_prod <- as.matrix(coef(cv_enet_prod, s = "lambda.min"))

eta_train <- as.numeric(Z_train %*% beta_enet_prod)

# Compute Breslow Cumulative Baseline Hazard on Dev
t_tr <- surv_train[, 1]; d_tr <- surv_train[, 2]
event_times <- sort(unique(t_tr[d_tr == 1]))
dH0 <- numeric(length(event_times))
for (j in seq_along(event_times)) {
  tj <- event_times[j]
  dj <- sum(t_tr == tj & d_tr == 1)
  risk_sum <- sum(exp(eta_train[t_tr >= tj]))
  dH0[j] <- dj / max(risk_sum, 1e-8)
}
H0 <- cumsum(dH0)
h0_step <- stepfun(event_times, c(0, H0))

# Save Frozen Production Model Bundle
frozen_bundle <- list(
  model_name = "ElasticNet-Cox (alpha = 0.5)",
  training_cohort_n = 285,
  training_deaths = sum(dev_events$os_event == 1),
  training_censored = sum(dev_events$os_event == 0),
  cv_fit = cv_enet_prod,
  lambda_min = lambda_star,
  coefficients = beta_enet_prod,
  mofa_fit = mofa_fit_tr,
  W_sub_list = W_sub_list,
  mu_tr_list = mu_tr_list,
  Omega_tr = Omega_tr,
  inv_omega = inv_omega,
  h0_step = h0_step,
  event_times = event_times,
  median_train_risk = median(eta_train),
  train_risk_mean = mean(eta_train),
  train_risk_sd = sd(eta_train),
  selected_features = list(RNA = rna_sel_genes, Methylation = meth_sel_probes, CNV = cnv_sel_genes, SNV = snv_sel_genes),
  combat_params_rna = combat_params_rna,
  combat_params_meth = combat_params_meth,
  tr_norm_factors_median = median(tr_norm_factors),
  tr_meth_meds = tr_meth_meds,
  annot_unique = annot_unique
)
saveRDS(frozen_bundle, "./FINAL_285_FROZEN_MODEL_BUNDLE.rds")
cat("Frozen production model bundle saved to FINAL_285_FROZEN_MODEL_BUNDLE.rds\n")

# 4. UNLOCK AND EVALUATE TRUE INDEPENDENT N = 65 HOLDOUT COHORT
cat("\n>>> UNLOCKING AND EVALUATING INDEPENDENT N = 65 HOLDOUT COHORT <<<\n")

# Project RNA Holdout
rna_raw_val <- rna_mat_all[, test_ids, drop = FALSE]
raw_val_mapped <- rna_raw_val[annot_unique$ORIG_ENS, , drop = FALSE]
rownames(raw_val_mapped) <- annot_unique$SYMBOL
lib_sizes_val <- colSums(rna_raw_val)
eff_lib_val <- lib_sizes_val * median(tr_norm_factors)
cpm_val_raw <- t(t((raw_val_mapped + 2) / (eff_lib_val + 2)) * 1e6)
log2_cpm_val <- log2(cpm_val_raw)
batch_val_rna <- substr(rna_aliqs[colnames(log2_cpm_val)], 22, 25)
rna_test_mat <- apply_combat_frozen(log2_cpm_val, batch_val_rna, combat_params_rna)[rna_sel_genes, , drop = FALSE]

# Project Methylation Holdout
meth_raw_val <- meth_mat_all[keep_probes_m, test_ids, drop = FALSE]
na_idx_val <- which(is.na(meth_raw_val), arr.ind = TRUE)
if (nrow(na_idx_val) > 0) meth_raw_val[na_idx_val] <- tr_meth_meds[na_idx_val[, 1]]
m_val <- beta2M(meth_raw_val)
m_val_sel <- m_val[meth_sel_probes, , drop = FALSE]
batch_val_meth <- substr(meth_aliqs[colnames(m_val_sel)], 22, 25)
meth_test_mat <- apply_combat_frozen(m_val_sel, batch_val_meth, combat_params_meth)

# CNV & SNV Holdout
cnv_test_mat <- cnv_full[cnv_sel_genes, test_ids, drop = FALSE]
snv_test_mat <- snv_full[snv_sel_genes, test_ids, drop = FALSE]

# Forward Project Holdout into MOFA Factor Space
U_test <- matrix(0, nrow = length(test_ids), ncol = 15)
for (vname in c("RNA", "Methylation", "CNV", "SNV")) {
  X_test_v <- switch(vname, RNA = rna_test_mat, Methylation = meth_test_mat, CNV = cnv_test_mat, SNV = snv_test_mat)
  W_sub    <- W_sub_list[[vname]]
  mu_tr_v  <- mu_tr_list[[vname]]
  
  X_test_cent <- X_test_v - mu_tr_v
  U_test <- U_test + t(X_test_cent) %*% W_sub
}

Z_test <- U_test %*% inv_omega
colnames(Z_test) <- paste0("Factor", 1:15)
rownames(Z_test) <- test_ids

# Compute Holdout Risk Scores
risk_holdout <- as.numeric(Z_test %*% beta_enet_prod)

# Compute Survival Probabilities
eval_times <- c(365, 1095, 1825)
surv_mat_holdout <- matrix(NA, nrow = length(test_ids), ncol = 3)
colnames(surv_mat_holdout) <- paste0("t_", eval_times)
for (k in seq_along(eval_times)) {
  t_star <- eval_times[k]
  surv_mat_holdout[, k] <- exp(-h0_step(t_star) * exp(risk_holdout))
}
surv_mat_holdout <- pmin(pmax(surv_mat_holdout, 0), 1)
for (k in 2:3) surv_mat_holdout[, k] <- pmin(surv_mat_holdout[, k], surv_mat_holdout[, k - 1])

# Save Holdout Predictions
holdout_preds_dt <- data.table(
  patient_id = test_ids,
  os_time = test_events$os_time,
  os_event = test_events$os_event,
  vital_status = test_events$vital_status,
  risk_score = risk_holdout,
  risk_group = ifelse(risk_holdout >= median(eta_train), "High", "Low"),
  surv_prob_1yr = surv_mat_holdout[, 1],
  surv_prob_3yr = surv_mat_holdout[, 2],
  surv_prob_5yr = surv_mat_holdout[, 3]
)
fwrite(holdout_preds_dt, "./FINAL_285_TRUE_HOLDOUT_PREDICTIONS.csv")
cat("Saved FINAL_285_TRUE_HOLDOUT_PREDICTIONS.csv\n")

# Calculate Holdout Metrics
c_holdout <- calc_cindex(surv_test, risk_holdout)

# 10,000 Bootstrap Confidence Interval
cat("Calculating 10,000 bootstrap iterations for C-index CI...\n")
set.seed(42003)
B <- 10000
boot_c <- numeric(B)
n_test <- length(test_ids)
for (b in 1:B) {
  b_idx <- sample(n_test, replace = TRUE)
  b_surv <- surv_test[b_idx]
  if (sum(b_surv[, 2]) > 0) {
    boot_c[b] <- calc_cindex(b_surv, risk_holdout[b_idx])
  } else {
    boot_c[b] <- NA
  }
}
boot_c <- boot_c[!is.na(boot_c)]
ci_lower <- as.numeric(quantile(boot_c, 0.025))
ci_upper <- as.numeric(quantile(boot_c, 0.975))
boot_se  <- sd(boot_c)

# Brier and IBS
brier_res <- calc_true_ipcw_brier(surv_train, surv_test, surv_mat_holdout)

# Kaplan-Meier Risk Group Analysis
km_fit <- survfit(surv_test ~ risk_group, data = holdout_preds_dt)
cox_km <- tryCatch(coxph(surv_test ~ risk_group, data = holdout_preds_dt), error = function(e) NULL)
km_diff <- survdiff(surv_test ~ risk_group, data = holdout_preds_dt)
p_logrank <- 1 - pchisq(km_diff$chisq, df = length(km_diff$n) - 1)
hr_km <- if (!is.null(cox_km)) exp(coef(cox_km)) else NA
hr_lower <- if (!is.null(cox_km)) exp(confint(cox_km))[1] else NA
hr_upper <- if (!is.null(cox_km)) exp(confint(cox_km))[2] else NA

high_n <- sum(holdout_preds_dt$risk_group == "High")
low_n  <- sum(holdout_preds_dt$risk_group == "Low")
high_events <- sum(holdout_preds_dt$risk_group == "High" & holdout_preds_dt$os_event == 1)
low_events  <- sum(holdout_preds_dt$risk_group == "Low"  & holdout_preds_dt$os_event == 1)

# Metrics Table
metrics_dt <- data.table(
  Metric = c(
    "Holdout_C_Index",
    "Holdout_C_Index_95CI_Lower",
    "Holdout_C_Index_95CI_Upper",
    "Holdout_C_Index_Bootstrap_SE",
    "Cohort_N",
    "Death_Events",
    "Censored_Patients",
    "Event_Rate_Pct",
    "Risk_Score_Mean",
    "Risk_Score_SD",
    "Risk_Score_Median",
    "Risk_Group_High_N",
    "Risk_Group_High_Events",
    "Risk_Group_Low_N",
    "Risk_Group_Low_Events",
    "Kaplan_Meier_Hazard_Ratio",
    "Kaplan_Meier_HR_95CI_Lower",
    "Kaplan_Meier_HR_95CI_Upper",
    "Kaplan_Meier_LogRank_PValue",
    "Surv_Prob_1y_Mean",
    "Surv_Prob_3y_Mean",
    "Surv_Prob_5y_Mean",
    "IPCW_Brier_Score_1y",
    "IPCW_Brier_Score_3y",
    "IPCW_Brier_Score_5y",
    "Integrated_Brier_Score_IBS"
  ),
  Value = c(
    c_holdout,
    ci_lower,
    ci_upper,
    boot_se,
    length(test_ids),
    sum(test_events$os_event == 1),
    sum(test_events$os_event == 0),
    mean(test_events$os_event == 1) * 100,
    mean(risk_holdout),
    sd(risk_holdout),
    median(risk_holdout),
    high_n,
    high_events,
    low_n,
    low_events,
    hr_km,
    hr_lower,
    hr_upper,
    p_logrank,
    mean(surv_mat_holdout[, 1]),
    mean(surv_mat_holdout[, 2]),
    mean(surv_mat_holdout[, 3]),
    brier_res$brier_scores[1],
    brier_res$brier_scores[2],
    brier_res$brier_scores[3],
    brier_res$ibs
  )
)

fwrite(metrics_dt, "./FINAL_285_TRUE_HOLDOUT_METRICS.csv")
cat("Saved FINAL_285_TRUE_HOLDOUT_METRICS.csv\n")

cat("\n================================================================================\n")
cat("TRUE INDEPENDENT HOLDOUT RESULTS (N = 65, 19 Deaths, 46 Censored):\n")
cat(sprintf("  Harrell's C-index: %.4f (95%% CI: [%.4f, %.4f], SE: %.4f)\n", c_holdout, ci_lower, ci_upper, boot_se))
cat(sprintf("  Integrated Brier Score (IBS): %.4f\n", brier_res$ibs))
cat(sprintf("  Brier 1-Year: %.4f | 3-Year: %.4f | 5-Year: %.4f\n",
            brier_res$brier_scores[1], brier_res$brier_scores[2], brier_res$brier_scores[3]))
cat(sprintf("  Kaplan-Meier HR: %.3f (95%% CI: [%.3f, %.3f], Log-Rank p = %.5f)\n",
            hr_km, hr_lower, hr_upper, p_logrank))
cat(sprintf("  Risk Groups: High N = %d (%d deaths) vs Low N = %d (%d deaths)\n",
            high_n, high_events, low_n, low_events))
cat("================================================================================\n")
