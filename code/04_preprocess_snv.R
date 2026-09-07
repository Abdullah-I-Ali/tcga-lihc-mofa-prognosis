# ==============================================================================
# NOTE: LEVEL 3 RAW RECONSTRUCTION PIPELINE
# This script filters somatic mutations to 9 functional classes and constructs
# the binary mutation matrix. It requires the raw input somatic MAF file:
# snv_mutation_raw.rds (acquired from NCI GDC; see docs/data_access.md).
# ==============================================================================
# ==============================================================================
# GATE 2B.4 — CLEAN SOMATIC MUTATION (SNV) PREPROCESSING PIPELINE
# ==============================================================================
# Modality: Somatic Single Nucleotide Variants & Small Indels (MAF format)
# Raw Input: snv_mutation_raw.rds
# Cohort: N = 350 Locked Primary Modeling Patients
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
})

data_dir <- ".""

cat("============================================================\n")
cat("  GATE 2B.4 CLEAN SNV PREPROCESSING PIPELINE\n")
cat("  Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n\n")

# ------------------------------------------------------------------------------
# 1. LOAD CANONICAL TRUTH TABLE & LOCKED COHORT
# ------------------------------------------------------------------------------
cat("1. Loading Canonical Cohort Truth Table...\n")
truth_df <- fread(file.path(data_dir, "GATE_2A6_FINAL_VERIFIED_PATIENT_TRUTH.csv"))
locked_pts <- sort(truth_df[PRIMARY_INCLUDED == TRUE, patient_id])
cat(sprintf("   - Locked Primary Modeling Cohort: N = %d Patients\n\n", length(locked_pts)))

# ------------------------------------------------------------------------------
# 2. LOAD RAW SNV MAF DATA & PARSE BARCODES
# ------------------------------------------------------------------------------
cat("2. Loading Raw SNV MAF Data and Parsing Sample Barcodes...\n")
snv_raw <- as.data.table(readRDS(file.path(data_dir, "snv_mutation_raw.rds")))
n_raw_muts <- nrow(snv_raw)
cat(sprintf("   - Raw Mutation Records: %d rows\n", n_raw_muts))

snv_raw[, Patient_ID := substr(Tumor_Sample_Barcode, 1, 12)]
snv_raw[, Sample_Type := substr(Tumor_Sample_Barcode, 14, 15)]

sample_type_names <- c(
  "01" = "Primary Solid Tumor",
  "02" = "Recurrent Solid Tumor",
  "10" = "Blood Derived Normal",
  "11" = "Solid Tissue Normal"
)
snv_raw[, Sample_Type_Name := sample_type_names[Sample_Type]]
snv_raw[is.na(Sample_Type_Name), Sample_Type_Name := "Other / Control"]

# ------------------------------------------------------------------------------
# 3. BUILD COMPLETE SNV SAMPLE MANIFEST
# ------------------------------------------------------------------------------
cat("3. Generating Complete SNV Sample Manifest...\n")

manifest_dt <- unique(snv_raw[, .(
  TUMOR_SAMPLE_BARCODE = Tumor_Sample_Barcode,
  PATIENT_ID = Patient_ID,
  SAMPLE_TYPE_CODE = Sample_Type,
  SAMPLE_TYPE_NAME = Sample_Type_Name
)])

manifest_dt[, IN_LOCKED_PRIMARY := PATIENT_ID %in% locked_pts]
manifest_dt[, SELECTION_STATUS := "EXCLUDED: Non-primary or non-locked cohort"]
manifest_dt[SAMPLE_TYPE_CODE == "01" & IN_LOCKED_PRIMARY == TRUE,
            SELECTION_STATUS := "SELECTED: Primary Tumor Aliquot for Locked Cohort"]

# Resolve multiple aliquots if any
aliq_counts <- manifest_dt[SELECTION_STATUS == "SELECTED: Primary Tumor Aliquot for Locked Cohort", .N, by = PATIENT_ID]
if (any(aliq_counts$N > 1)) {
  # deterministic deduplication
  dup_pts <- aliq_counts[N > 1, PATIENT_ID]
  for (dp in dup_pts) {
    aliqs <- sort(manifest_dt[PATIENT_ID == dp & SELECTION_STATUS == "SELECTED: Primary Tumor Aliquot for Locked Cohort", TUMOR_SAMPLE_BARCODE])
    manifest_dt[PATIENT_ID == dp & TUMOR_SAMPLE_BARCODE != aliqs[1], SELECTION_STATUS := "EXCLUDED: Duplicate primary aliquot"]
  }
}

write.csv(manifest_dt, file.path(data_dir, "GATE_2B4_SNV_SAMPLE_MANIFEST.csv"), row.names = FALSE)
cat("   ✓ Saved: GATE_2B4_SNV_SAMPLE_MANIFEST.csv\n\n")

# Filter SNV records to selected primary tumor barcodes
selected_barcodes <- manifest_dt[SELECTION_STATUS == "SELECTED: Primary Tumor Aliquot for Locked Cohort", TUMOR_SAMPLE_BARCODE]
snv_selected <- snv_raw[Tumor_Sample_Barcode %in% selected_barcodes]
cat(sprintf("   - Retained Primary Tumor Mutation Records: %d rows across %d Unique Patients (1:1 Mapped)\n\n",
            nrow(snv_selected), length(unique(snv_selected$Patient_ID))))

# ------------------------------------------------------------------------------
# 4. MUTATION-LEVEL FUNCTIONAL CLASSIFICATION & QC
# ------------------------------------------------------------------------------
cat("4. Performing Mutation-Level Functional Classification...\n")

functional_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
  "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation", "Translation_Start_Site"
)

snv_selected[, is_functional := Variant_Classification %in% functional_classes]

n_func <- sum(snv_selected$is_functional)
n_nonfunc <- sum(!snv_selected$is_functional)

cat(sprintf("   - Functional Somatic Mutations Retained: %d (%.2f%%)\n", n_func, (n_func / nrow(snv_selected)) * 100))
cat(sprintf("   - Silent / Non-Functional Excluded:     %d (%.2f%%)\n\n", n_nonfunc, (n_nonfunc / nrow(snv_selected)) * 100))

# Filter to functional mutations
snv_functional <- snv_selected[is_functional == TRUE]

# ------------------------------------------------------------------------------
# 5. CONSTRUCT BINARY MUTATION MATRIX & SPARSITY METRICS
# ------------------------------------------------------------------------------
cat("5. Constructing Binary Somatic Mutation Matrices...\n")

# Aggregate binary indicator per gene x patient
pt_gene_pairs <- unique(snv_functional[, .(Hugo_Symbol, Patient_ID)])
pt_gene_pairs[, mutated := 1L]

# Dense matrix casting
snv_cast <- dcast(pt_gene_pairs, Hugo_Symbol ~ Patient_ID, value.var = "mutated", fill = 0L)
gene_symbols <- snv_cast$Hugo_Symbol
snv_cast[, Hugo_Symbol := NULL]
full_snv_matrix <- as.matrix(snv_cast)
rownames(full_snv_matrix) <- gene_symbols

# Ensure all 350 patients are present and sorted
missing_pts <- setdiff(locked_pts, colnames(full_snv_matrix))
if (length(missing_pts) > 0) {
  for (mp in missing_pts) {
    full_snv_matrix <- cbind(full_snv_matrix, setNames(matrix(0L, nrow = nrow(full_snv_matrix), ncol = 1), mp))
  }
}
full_snv_matrix <- full_snv_matrix[, locked_pts, drop = FALSE]

# Sparsity & prevalence metrics
total_elements <- length(full_snv_matrix)
zero_elements  <- sum(full_snv_matrix == 0L)
overall_sparsity <- (zero_elements / total_elements) * 100

cat(sprintf("   - Full Binary SNV Matrix: %d Genes x %d Patients\n", nrow(full_snv_matrix), ncol(full_snv_matrix)))
cat(sprintf("   - Overall Sparsity: %.4f%% zeros\n\n", overall_sparsity))

# ------------------------------------------------------------------------------
# 6. FEATURE SELECTION FOR MOFA2 (TOP 500 RECURRENT MUTATED GENES)
# ------------------------------------------------------------------------------
cat("6. Selecting Top 500 Recurrent Mutated Genes for MOFA2 (Bernoulli Likelihood)...\n")

gene_recurrence <- rowSums(full_snv_matrix)
top_500_idx <- order(gene_recurrence, decreasing = TRUE)[1:500]
top_500_genes <- sort(rownames(full_snv_matrix)[top_500_idx]) # Deterministic alphabetical sort
mofa_snv_matrix <- full_snv_matrix[top_500_genes, , drop = FALSE]

top_sparsity <- (sum(mofa_snv_matrix == 0L) / length(mofa_snv_matrix)) * 100
cat(sprintf("   - MOFA Binary SNV Matrix: %d Genes x %d Patients\n", nrow(mofa_snv_matrix), ncol(mofa_snv_matrix)))
cat(sprintf("   - Top 500 Recurrent Gene Sparsity: %.4f%% zeros\n\n", top_sparsity))

# ------------------------------------------------------------------------------
# 7. SERIALIZE DELIVERABLE ARTIFACT
# ------------------------------------------------------------------------------
cat("7. Serializing GATE_2B4_SNV_MATRIX.rds...\n")
snv_output_obj <- list(
  mofa_features = mofa_snv_matrix,
  full_matrix = full_snv_matrix,
  selected_genes = top_500_genes,
  mutation_summary = snv_selected[, .N, by = .(Variant_Classification, is_functional)],
  patient_tmb = snv_functional[, .(functional_mutations = .N), by = Patient_ID],
  metadata = list(
    modality = "Somatic Single Nucleotide Variants & Indels (Binary)",
    data_format = "MAF",
    functional_classes = functional_classes,
    top_k_features = 500,
    likelihood_model = "Bernoulli",
    timestamp = Sys.time()
  )
)

saveRDS(snv_output_obj, file.path(data_dir, "GATE_2B4_SNV_MATRIX.rds"))
cat("   ✓ Saved: GATE_2B4_SNV_MATRIX.rds\n\n")

# ------------------------------------------------------------------------------
# 8. INDEPENDENT VERIFICATION ASSERTIONS
# ------------------------------------------------------------------------------
cat("8. Running Automated Verification Assertions...\n")

assertions <- list()

assertions[[1]] <- list(
  CHECK_ID = "EXACT_N_350_PATIENTS",
  RESULT = if (ncol(mofa_snv_matrix) == 350 && ncol(full_snv_matrix) == 350) "PASS" else "FAIL",
  MESSAGE = sprintf("Expected 350, got full=%d, mofa=%d", ncol(full_snv_matrix), ncol(mofa_snv_matrix))
)

assertions[[2]] <- list(
  CHECK_ID = "ZERO_DUPLICATE_PATIENTS",
  RESULT = if (length(unique(colnames(mofa_snv_matrix))) == 350) "PASS" else "FAIL",
  MESSAGE = "Duplicate patient barcodes detected"
)

assertions[[3]] <- list(
  CHECK_ID = "PATIENTS_MATCH_LOCKED_COHORT",
  RESULT = if (identical(colnames(mofa_snv_matrix), locked_pts) && identical(colnames(full_snv_matrix), locked_pts)) "PASS" else "FAIL",
  MESSAGE = "Patient IDs mismatch canonical locked cohort"
)

assertions[[4]] <- list(
  CHECK_ID = "EXACT_TOP_GENES_500",
  RESULT = if (nrow(mofa_snv_matrix) == 500) "PASS" else "FAIL",
  MESSAGE = sprintf("Expected 500 genes, got %d", nrow(mofa_snv_matrix))
)

assertions[[5]] <- list(
  CHECK_ID = "BINARY_VALUES_ONLY",
  RESULT = if (all(mofa_snv_matrix %in% c(0L, 1L)) && all(full_snv_matrix %in% c(0L, 1L))) "PASS" else "FAIL",
  MESSAGE = "Non-binary values present in SNV matrix"
)

assertions[[6]] <- list(
  CHECK_ID = "ZERO_NAN_INF_NA",
  RESULT = if (sum(is.na(mofa_snv_matrix) | is.nan(mofa_snv_matrix) | is.infinite(mofa_snv_matrix)) == 0) "PASS" else "FAIL",
  MESSAGE = "Non-finite values present in matrix"
)

assertions[[7]] <- list(
  CHECK_ID = "REPRODUCIBLE_FEATURE_ORDER",
  RESULT = if (identical(rownames(mofa_snv_matrix), sort(rownames(mofa_snv_matrix)))) "PASS" else "FAIL",
  MESSAGE = "Features not deterministically ordered"
)

excluded_hypermutators <- truth_df[hypermutator == TRUE, patient_id]
assertions[[8]] <- list(
  CHECK_ID = "HYPERMUTATORS_EXCLUDED",
  RESULT = if (length(intersect(colnames(mofa_snv_matrix), excluded_hypermutators)) == 0) "PASS" else "FAIL",
  MESSAGE = "Hypermutators present in primary cohort"
)

assert_df <- rbindlist(assertions)
write.csv(assert_df, file.path(data_dir, "GATE_2B4_SNV_VERIFICATION.csv"), row.names = FALSE)
cat("   ✓ Saved: GATE_2B4_SNV_VERIFICATION.csv\n\n")
print(assert_df)

cat("\n============================================================\n")
cat("  GATE 2B.4 CLEAN SNV PREPROCESSING COMPLETE\n")
cat("  Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")
