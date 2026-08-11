#!/usr/bin/env Rscript
# Diagnose the "unclassified" bucket in report 40 (TCGA classification).
# Run from the ATTEND repo root, WHERE THE DATA EXISTS (HPC, or after syncing
# data/gistic/ + output/clean_data/ + the gianlu clinical file down locally):
#   Rscript code/diagnose_unclassified.R
# It reproduces report 40's objects up through the cascade and prints where
# patients are lost between the master, the GISTIC matrix, and the crosswalk.

suppressPackageStartupMessages({
  library(tidyverse); library(here)
})

source(here("code", "attend_harmonise.R"))
source(here("code", "attend_classes.R"))
source(here("code", "attend_io.R"))
source(here("code", "load_clinical.R"))
source(here("code", "load_wes_results.R"))

k_cnv <- 4  # doesn't affect the unclassified set; only how clustered samples split

# --- rebuild the report-09 objects -----------------------------------------
master <- tryCatch(read_intermediate("attend_master_joined") |> add_molecular_classes(),
                   error = function(e) { message("master unavailable: ", conditionMessage(e)); NULL })
cw     <- tryCatch(build_barcode_pid(load_gianlu_clinical_data(), make_id_cfg()),
                   error = function(e) { message("crosswalk unavailable: ", conditionMessage(e)); NULL })
cnv_mat <- tryCatch(load_gistic_thresholded(), error = function(e) NULL)

stopifnot("master missing"   = !is.null(master),
          "crosswalk missing" = !is.null(cw),
          "GISTIC matrix missing" = !is.null(cnv_mat))

cl <- cluster_arm_matrix(cnv_mat, k = k_cnv, method = "complete", dist_method = "euclidean")
high_cluster <- cnv_high_cluster_burden(cl)

cnv_pid <- cl$assignment |>
  inner_join(cw, by = c("ID" = "barcode")) |>
  count(pid, cnv_cluster) |>
  group_by(pid) |> slice_max(n, n = 1, with_ties = FALSE) |> ungroup() |>
  transmute(pid, cnv_high = cnv_cluster == high_cluster)

classed <- master |> left_join(cnv_pid, by = "pid") |> add_tcga_class()

# --- the shrink points ------------------------------------------------------
cat("\n================ WHERE SAMPLES ARE LOST ================\n")
cat(sprintf("master patients (pid)            : %d\n", nrow(master)))
cat(sprintf("samples in GISTIC matrix         : %d\n", nrow(cl$assignment)))
gistic_ids <- unique(cl$assignment$ID)
matched    <- intersect(gistic_ids, cw$barcode)
cat(sprintf("GISTIC IDs matching crosswalk    : %d of %d   <-- key number\n",
            length(matched), length(gistic_ids)))
cat(sprintf("patients with a CNV cluster      : %d\n", sum(!is.na(classed$cnv_high))))

cat("\n---- integrated class breakdown (NA = unclassified) ----\n")
print(table(classed$tcga_class, useNA = "always"))

# MMRd status among the unclassified — confirms they're the non-MMRd, non-CNV set
uncl <- classed |> filter(is.na(tcga_class))
cat(sprintf("\nunclassified patients            : %d\n", nrow(uncl)))
cat(sprintf("  of which MMR/MSI status known  : %d\n",
            sum(!is.na(uncl$MSI_class) | !is.na(uncl$MMR_class))))
cat(sprintf("  of which have cnv_high (should be 0): %d\n", sum(!is.na(uncl$cnv_high))))

# --- if the key number is small, eyeball the ID-format mismatch -------------
if (length(matched) < 0.9 * length(gistic_ids)) {
  cat("\n!!! ID MISMATCH suspected — GISTIC IDs the crosswalk does NOT contain:\n")
  print(utils::head(setdiff(gistic_ids, cw$barcode), 10))
  cat("... vs what the crosswalk (TUMOR_BARCODE) expects:\n")
  print(utils::head(cw$barcode, 10))
  cat("\nFix: adjust attend_cnv$seg$id_strip in code/attend_classes.R so the\n",
      "stripped .seg/GISTIC name equals TUMOR_BARCODE, then re-run.\n", sep = "")
}
cat("\n========================================================\n")
