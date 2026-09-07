suppressPackageStartupMessages({
  library(data.table)
  library(caret)
})

dev_ids_dt <- fread("./TRUE_DEV_285_IDS.csv")
dev_ids <- sort(dev_ids_dt$patient_id)
stopifnot(length(dev_ids) == 285)

event_dt <- fread("./TRUE_DEV_285_EVENT_VECTOR.csv")
event_dt <- event_dt[match(dev_ids, patient_id)]
stopifnot(identical(event_dt$patient_id, dev_ids))
stopifnot(sum(event_dt$os_event == 1) == 103)
stopifnot(sum(event_dt$os_event == 0) == 182)

# Deterministic 5 repeats x 5 folds = 25 outer runs
set.seed(42)
multi_folds <- createMultiFolds(as.factor(event_dt$os_event), k = 5, times = 5)
stopifnot(length(multi_folds) == 25)

fold_names <- names(multi_folds)
cat("Generated 25 stratified outer folds on N = 285.\n")

records <- list()
for (i in seq_along(multi_folds)) {
  fname <- fold_names[i]
  parts <- strsplit(fname, "\\.")[[1]]
  fold_num <- as.integer(sub("Fold", "", parts[1]))
  rep_num  <- as.integer(sub("Rep", "", parts[2]))
  
  tr_idx  <- multi_folds[[i]]
  val_idx <- setdiff(1:285, tr_idx)
  
  tr_ids  <- dev_ids[tr_idx]
  val_ids <- dev_ids[val_idx]
  
  stopifnot(length(intersect(tr_ids, val_ids)) == 0)
  stopifnot(length(union(tr_ids, val_ids)) == 285)
  stopifnot(length(val_ids) >= 56 && length(val_ids) <= 58)
  stopifnot(length(tr_ids) + length(val_ids) == 285)
  
  tr_events  <- event_dt[patient_id %in% tr_ids, os_event]
  val_events <- event_dt[patient_id %in% val_ids, os_event]
  
  records[[i]] <- data.table(
    run_index = i,
    repeat_id = rep_num,
    fold_id = fold_num,
    run_id = sprintf("Fold%d.Rep%02d", fold_num, rep_num),
    outer_seed = 42,
    mofa_seed = 100000 + i,
    inner_cv_seed = 200000 + i,
    lasso_seed = 300000 + i,
    xgboost_seed = 400000 + i,
    rsf_seed = 500000 + i,
    bootstrap_seed = 600000 + i,
    train_N = length(tr_ids),
    validation_N = length(val_ids),
    train_deaths = sum(tr_events == 1),
    train_censored = sum(tr_events == 0),
    val_deaths = sum(val_events == 1),
    val_censored = sum(val_events == 0),
    disjoint_pass = TRUE,
    complete_pass = TRUE
  )
}

fold_matrix <- rbindlist(records)
fwrite(fold_matrix, "./TRUE_DEV_285_SEED_FOLD_MATRIX.csv")
cat("Saved ./TRUE_DEV_285_SEED_FOLD_MATRIX.csv\n")
print(head(fold_matrix, 5))
