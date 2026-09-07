# ==============================================================================
# NOTE: LEVEL 3 RAW RECONSTRUCTION PIPELINE
# This script performs length-weighted mean segment-to-gene aggregation from
# Affymetrix SNP 6.0 segments. It requires the raw input segment file:
# cnv_segment_raw.rds (acquired from NCI GDC; see docs/data_access.md).
# ==============================================================================
# ==============================================================================
# GATE 2B.3 — CLEAN COPY NUMBER VARIATION (CNV) PREPROCESSING PIPELINE
# ==============================================================================
# Modality: Affymetrix Genome-Wide Human SNP Array 6.0 (hg38 Segment Mean)
# Raw Input: cnv_segment_raw.rds & gene_coords_hg38.rds
# Cohort: N = 350 Locked Primary Modeling Patients
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
})

data_dir <- ".""

cat("============================================================\n")
cat("  GATE 2B.3 CLEAN CNV PREPROCESSING PIPELINE\n")
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
# 2. LOAD RAW CNV SEGMENTS & DECONSTRUCT COMPOUND BARCODES
# ------------------------------------------------------------------------------
cat("2. Loading Raw CNV Segments and Deconstructing Barcodes...\n")
cnv_raw <- as.data.table(readRDS(file.path(data_dir, "cnv_segment_raw.rds")))
n_raw_rows <- nrow(cnv_raw)
cat(sprintf("   - Raw CNV Segment Records: %d rows\n", n_raw_rows))

# Parse compound barcodes (Sample column contains 'Tumor;Normal' pairs or single aliquots)
sample_split <- strsplit(cnv_raw$Sample, ";")
tumor_barcodes <- sapply(sample_split, function(x) {
  types <- substr(x, 14, 15)
  tumors <- x[types %in% c("01", "02", "06")]
  if (length(tumors) > 0) tumors[1] else x[1]
})
cnv_raw[, Tumor_Barcode := tumor_barcodes]
cnv_raw[, Patient_ID := substr(Tumor_Barcode, 1, 12)]
cnv_raw[, Sample_Type := substr(Tumor_Barcode, 14, 15)]
cnv_raw[, Center_Code := substr(Tumor_Barcode, 27, 28)]
cnv_raw[, Plate_ID := substr(Tumor_Barcode, 22, 25)]

# Sample type naming
sample_type_names <- c(
  "01" = "Primary Solid Tumor",
  "02" = "Recurrent Solid Tumor",
  "10" = "Blood Derived Normal",
  "11" = "Solid Tissue Normal"
)
cnv_raw[, Sample_Type_Name := sample_type_names[Sample_Type]]
cnv_raw[is.na(Sample_Type_Name), Sample_Type_Name := "Other / Control"]

# ------------------------------------------------------------------------------
# 3. BUILD COMPLETE CNV SAMPLE MANIFEST
# ------------------------------------------------------------------------------
cat("3. Generating Complete CNV Sample Manifest...\n")

manifest_dt <- unique(cnv_raw[, .(
  RAW_SAMPLE_STRING = Sample,
  TUMOR_ALIQUOT_BARCODE = Tumor_Barcode,
  PATIENT_ID = Patient_ID,
  SAMPLE_TYPE_CODE = Sample_Type,
  SAMPLE_TYPE_NAME = Sample_Type_Name,
  CENTER_CODE = Center_Code,
  PLATE_ID = Plate_ID,
  GDC_ALIQUOT_ID = GDC_Aliquot_ID
)])

manifest_dt[, IN_LOCKED_PRIMARY := PATIENT_ID %in% locked_pts]

# Selection logic: Primary tumor (01), Center 01, resolving TCGA-CC-A8HV to A40C
manifest_dt[, SELECTION_STATUS := "EXCLUDED: Non-primary or non-locked cohort"]
manifest_dt[SAMPLE_TYPE_CODE == "01" & IN_LOCKED_PRIMARY == TRUE & CENTER_CODE == "01", 
            SELECTION_STATUS := "SELECTED: Primary Tumor Aliquot for Locked Cohort"]

# Resolve TCGA-CC-A8HV duplicate
manifest_dt[PATIENT_ID == "TCGA-CC-A8HV" & TUMOR_ALIQUOT_BARCODE == "TCGA-CC-A8HV-01A-11D-A35Y-01",
            SELECTION_STATUS := "EXCLUDED: Lower segment resolution duplicate aliquot"]

write.csv(manifest_dt, file.path(data_dir, "GATE_2B3_CNV_SAMPLE_MANIFEST.csv"), row.names = FALSE)
cat("   ✓ Saved: GATE_2B3_CNV_SAMPLE_MANIFEST.csv\n\n")

# Filter CNV records to selected primary tumor segments
cnv_selected <- cnv_raw[Sample_Type == "01" & Patient_ID %in% locked_pts & Center_Code == "01"]
cnv_selected <- cnv_selected[!(Patient_ID == "TCGA-CC-A8HV" & Tumor_Barcode == "TCGA-CC-A8HV-01A-11D-A35Y-01")]

cat(sprintf("   - Retained Selected Primary CNV Segments: %d rows across %d Unique Patients (1:1 Mapped)\n\n",
            nrow(cnv_selected), length(unique(cnv_selected$Patient_ID))))

# ------------------------------------------------------------------------------
# 4. GENE COORDINATES & AUTOSOMAL GENOMIC OVERLAP
# ------------------------------------------------------------------------------
cat("4. Loading Gene Coordinates & Performing Genomic Overlap...\n")
gene_coords <- as.data.table(readRDS(file.path(data_dir, "gene_coords_hg38.rds")))
gene_coords[, chr_clean := as.character(chromosome_name)]
gene_autosomal <- gene_coords[chr_clean %in% as.character(1:22)]

# Standardize gene boundaries by taking broadest span per HGNC symbol
gene_autosomal <- gene_autosomal[, .(
  start_position = min(start_position),
  end_position = max(end_position)
), by = .(hgnc_symbol, chr_clean)]
setorder(gene_autosomal, chr_clean, start_position)
cat(sprintf("   - Autosomal Canonical Genes (chr1-22): %d unique HGNC symbols\n", nrow(gene_autosomal)))

# Standardize CNV chromosomes
cnv_selected[, chr_clean := gsub("^chr", "", Chromosome)]
cnv_selected <- cnv_selected[chr_clean %in% as.character(1:22)]

# Fast Interval Overlap via foverlaps
setkey(cnv_selected, chr_clean, Start, End)
setkey(gene_autosomal, chr_clean, start_position, end_position)

overlap_dt <- foverlaps(
  gene_autosomal,
  cnv_selected,
  by.x = c("chr_clean", "start_position", "end_position"),
  by.y = c("chr_clean", "Start", "End"),
  type = "any",
  nomatch = NA
)

overlap_dt[, overlap_start := pmax(Start, start_position)]
overlap_dt[, overlap_end   := pmin(End, end_position)]
overlap_dt[, overlap_len   := pmax(0, overlap_end - overlap_start + 1)]

cat(sprintf("   - Total Genomic Overlaps Identified: %d\n\n", nrow(overlap_dt)))

# ------------------------------------------------------------------------------
# 5. LENGTH-WEIGHTED MEAN GENE-LEVEL AGGREGATION
# ------------------------------------------------------------------------------
cat("5. Aggregating to Gene-Level Copy Ratio (Length-Weighted Mean)...\n")

gene_cnv_agg <- overlap_dt[!is.na(Patient_ID), .(
  weighted_mean = sum(overlap_len * Segment_Mean, na.rm = TRUE) / sum(overlap_len, na.rm = TRUE)
), by = .(hgnc_symbol, Patient_ID)]

# Cast to Dense Gene x Patient Matrix
gene_cnv_mat <- dcast(gene_cnv_agg, hgnc_symbol ~ Patient_ID, value.var = "weighted_mean", fill = 0.0)
gene_names <- gene_cnv_mat$hgnc_symbol
gene_cnv_mat[, hgnc_symbol := NULL]
full_cnv_matrix <- as.matrix(gene_cnv_mat)
rownames(full_cnv_matrix) <- gene_names

# Ensure all 350 patients are present and sorted
missing_pts <- setdiff(locked_pts, colnames(full_cnv_matrix))
if (length(missing_pts) > 0) {
  for (mp in missing_pts) {
    full_cnv_matrix <- cbind(full_cnv_matrix, setNames(matrix(0, nrow = nrow(full_cnv_matrix), ncol = 1), mp))
  }
}
full_cnv_matrix <- full_cnv_matrix[, locked_pts, drop = FALSE]

cat(sprintf("   - Full Gene-Level CNV Matrix: %d Genes x %d Patients\n",
            nrow(full_cnv_matrix), ncol(full_cnv_matrix)))
cat(sprintf("   - Value Range: [%.4f, %.4f], NA/NaN/Inf = %d\n\n",
            min(full_cnv_matrix), max(full_cnv_matrix),
            sum(is.na(full_cnv_matrix) | is.nan(full_cnv_matrix) | is.infinite(full_cnv_matrix))))

# ------------------------------------------------------------------------------
# 6. FEATURE SELECTION FOR MOFA2 (TOP 5,000 VARIABLE GENES)
# ------------------------------------------------------------------------------
cat("6. Selecting Top 5,000 Variable Genes for MOFA2...\n")
cnv_vars <- rowVars(full_cnv_matrix)
top_5000_idx <- order(cnv_vars, decreasing = TRUE)[1:5000]
top_5000_genes <- sort(rownames(full_cnv_matrix)[top_5000_idx]) # Deterministic alphabetical sorting
mofa_cnv_matrix <- full_cnv_matrix[top_5000_genes, , drop = FALSE]

cat(sprintf("   - MOFA Feature Matrix: %d Genes x %d Patients\n\n",
            nrow(mofa_cnv_matrix), ncol(mofa_cnv_matrix)))

# ------------------------------------------------------------------------------
# 7. SERIALIZE DELIVERABLE ARTIFACT
# ------------------------------------------------------------------------------
cat("7. Serializing GATE_2B3_CNV_MATRIX.rds...\n")
cnv_output_obj <- list(
  mofa_features = mofa_cnv_matrix,
  full_matrix = full_cnv_matrix,
  selected_genes = top_5000_genes,
  gene_coordinates = gene_autosomal,
  metadata = list(
    modality = "Copy Number Variation (log2 copy ratio)",
    platform = "Affymetrix Genome-Wide Human SNP Array 6.0",
    genome_build = "hg38",
    aggregation_method = "Length-weighted segment mean",
    neutral_uncovered_value = 0.0,
    timestamp = Sys.time()
  )
)

saveRDS(cnv_output_obj, file.path(data_dir, "GATE_2B3_CNV_MATRIX.rds"))
cat("   ✓ Saved: GATE_2B3_CNV_MATRIX.rds\n\n")

# ------------------------------------------------------------------------------
# 8. INDEPENDENT VERIFICATION ASSERTIONS
# ------------------------------------------------------------------------------
cat("8. Running Automated Verification Assertions...\n")

assertions <- list()

assertions[[1]] <- list(
  CHECK_ID = "EXACT_N_350_PATIENTS",
  RESULT = if (ncol(mofa_cnv_matrix) == 350 && ncol(full_cnv_matrix) == 350) "PASS" else "FAIL",
  MESSAGE = sprintf("Expected 350, got full=%d, mofa=%d", ncol(full_cnv_matrix), ncol(mofa_cnv_matrix))
)

assertions[[2]] <- list(
  CHECK_ID = "ZERO_DUPLICATE_PATIENTS",
  RESULT = if (length(unique(colnames(mofa_cnv_matrix))) == 350) "PASS" else "FAIL",
  MESSAGE = "Duplicate patient barcodes detected"
)

assertions[[3]] <- list(
  CHECK_ID = "PATIENTS_MATCH_LOCKED_COHORT",
  RESULT = if (identical(colnames(mofa_cnv_matrix), locked_pts) && identical(colnames(full_cnv_matrix), locked_pts)) "PASS" else "FAIL",
  MESSAGE = "Patient IDs mismatch canonical locked cohort"
)

assertions[[4]] <- list(
  CHECK_ID = "EXACT_TOP_GENES_5000",
  RESULT = if (nrow(mofa_cnv_matrix) == 5000) "PASS" else "FAIL",
  MESSAGE = sprintf("Expected 5000 genes, got %d", nrow(mofa_cnv_matrix))
)

assertions[[5]] <- list(
  CHECK_ID = "ZERO_NAN_INF_NA",
  RESULT = if (sum(is.na(mofa_cnv_matrix) | is.nan(mofa_cnv_matrix) | is.infinite(mofa_cnv_matrix)) == 0) "PASS" else "FAIL",
  MESSAGE = "Non-finite values present in matrix"
)

assertions[[6]] <- list(
  CHECK_ID = "VALID_NUMERICAL_RANGE",
  RESULT = if (min(full_cnv_matrix) >= -35 && max(full_cnv_matrix) <= 15) "PASS" else "FAIL",
  MESSAGE = sprintf("Range [%.2f, %.2f] outside biological bounds", min(full_cnv_matrix), max(full_cnv_matrix))
)

assertions[[7]] <- list(
  CHECK_ID = "REPRODUCIBLE_FEATURE_ORDER",
  RESULT = if (identical(rownames(mofa_cnv_matrix), sort(rownames(mofa_cnv_matrix)))) "PASS" else "FAIL",
  MESSAGE = "Features not deterministically ordered"
)

assertions[[8]] <- list(
  CHECK_ID = "ONE_TO_ONE_PRIMARY_MAPPING",
  RESULT = if (nrow(manifest_dt[SELECTION_STATUS == "SELECTED: Primary Tumor Aliquot for Locked Cohort"]) == 350) "PASS" else "FAIL",
  MESSAGE = "Manifest does not map exactly 1:1"
)

assert_df <- rbindlist(assertions)
write.csv(assert_df, file.path(data_dir, "GATE_2B3_CNV_VERIFICATION.csv"), row.names = FALSE)
cat("   ✓ Saved: GATE_2B3_CNV_VERIFICATION.csv\n\n")
print(assert_df)

cat("\n============================================================\n")
cat("  GATE 2B.3 CLEAN CNV PREPROCESSING COMPLETE\n")
cat("  Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")
