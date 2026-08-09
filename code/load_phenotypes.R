# =============================================================================
# load_phenotypes.R  —  IHC / PHENOTYPE LOADER
#
# Sourced by reports 1 and 4. Reads the FlowPath per-cell IHC export (one
# CSV per patient under a directory tree) and reduces each patient to a single
# row of cell counts, keyed by `image_id` (mapped to pid via the image crosswalk).
# Patients whose folder doesn't hold exactly one CSV are skipped (possibly()).
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(fs)
  library(tidyverse)
  library(data.table)
})

# Cluster location of the per-patient FlowPath CSVs (override to run elsewhere).
FLOWPATH_DIR <- "/hpcnfs/techunits/imaging/work/ATTEND/FlowPath_csv/pheno_micro"

process_patient <- function(dir) {
  csv_path <- dir_ls(dir, glob = "*.csv")
  if (length(csv_path) != 1)
    stop(sprintf("found %d csv files (expected 1)", length(csv_path)))  # -> skipped by possibly()

  patient_id <- path_file(path_dir(csv_path))
  message(sprintf("LOADING PATIENT: %s", patient_id))

  ihc <- fread(csv_path) |>
    as_tibble() |>
    mutate(phenotype_clean = str_extract(phenotype, "(?<=\\().*?(?=\\))"))

  inside <- ihc |> filter(Out_of_annotation == FALSE)   # cells within the annotated region

  # PD-1 / PD-L1 staining is only present in some batches.
  has_pd1 <- any(str_detect(colnames(inside), "PD1"))
  # CD3 is part of the same multiplex panel but, like PD-1, absent in some batches.
  has_cd3 <- any(str_detect(colnames(inside), "CD3"))

  tibble(
    image_id              = patient_id,
    n_cells               = nrow(ihc),
    n_inside              = nrow(inside),
    n_tumor_total         = sum(str_detect(ihc$phenotype, "Tumor")),
    n_tumor_inside        = sum(str_detect(inside$phenotype, "Tumor")),
    n_cd45_inside         = sum(inside$CD45_sign == "+"),
    # CD3+ leukocytes and CD3+CD45+ double-positive T cells inside the annotation.
    n_cd3_inside          = if (has_cd3) sum(inside$CD3_sign == "+") else NA,
    n_cd3cd45_inside      = if (has_cd3) sum(inside$CD3_sign == "+" & inside$CD45_sign == "+") else NA,
    n_tumor_aridia_inside = sum(inside$ARID1A_sign == "+" &
                                  str_detect(inside$phenotype_clean, "Tumor")),
    n_tumor_pdl1_inside   = if (has_pd1) sum(inside$PDL1_sign == "+" &
                                               str_detect(inside$phenotype_clean, "Tumor")) else NA,
    n_cd45_pd1_inside     = if (has_pd1) sum(inside$PD1_sign  == "+" & inside$CD45_sign == "+") else NA,
    n_cd45_pdl1_inside    = if (has_pd1) sum(inside$PDL1_sign == "+" & inside$CD45_sign == "+") else NA
  )
}

load_ihc_data <- function(flowpath_dir = FLOWPATH_DIR) {
  safe_process <- possibly(process_patient, otherwise = NULL, quiet = FALSE)
  dir_ls(flowpath_dir) |>
    map(safe_process) |>
    list_rbind()
}

# --- Per-cell-type composition (inside the annotation) ----------------------
# For report 4: counts every CELL TYPE inside the annotation, per image, plus the
# three denominators (all cells / tumour cells / CD45+ cells inside) so downstream
# can normalise. "Cell type" = phenotype_clean (the parenthetical label in
# `phenotype`); switch to `phenotype` for the full label. Returns long rows:
# image_id, cell_type, n_cell, n_inside, n_tumor_inside, n_cd45_inside.
process_patient_celltypes <- function(dir) {
  csv_path <- dir_ls(dir, glob = "*.csv")
  if (length(csv_path) != 1)
    stop(sprintf("found %d csv files (expected 1)", length(csv_path)))
  image_id <- path_file(path_dir(csv_path))

  inside <- fread(csv_path) |>
    as_tibble() |>
    mutate(phenotype_clean = str_extract(phenotype, "(?<=\\().*?(?=\\))")) |>
    filter(Out_of_annotation == FALSE)

  has_cd3 <- any(str_detect(colnames(inside), "CD3"))

  inside |>
    count(cell_type = phenotype_clean, name = "n_cell") |>
    filter(!is.na(cell_type)) |>
    mutate(image_id         = image_id,
           n_inside         = nrow(inside),
           n_tumor_inside   = sum(str_detect(inside$phenotype, "Tumor")),
           n_cd45_inside    = sum(inside$CD45_sign == "+"),
           # CD3+CD45+ double-positive T cells inside — denominator for the immune metrics.
           n_cd3cd45_inside = if (has_cd3) sum(inside$CD3_sign == "+" & inside$CD45_sign == "+") else NA)
}

load_ihc_celltypes <- function(flowpath_dir = FLOWPATH_DIR) {
  safe_ct <- possibly(process_patient_celltypes, otherwise = NULL, quiet = FALSE)
  dir_ls(flowpath_dir) |>
    map(safe_ct) |>
    list_rbind()
}
