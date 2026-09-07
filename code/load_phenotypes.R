# =============================================================================
# load_phenotypes.R  —  IHC / PHENOTYPE LOADER
#
# Sourced by report 01, 20 and 32. Reads the FlowPath per-cell IHC export (one
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

# =============================================================================
# CACHING — why these loaders were slow, and what is done about it
#
# Both loaders walk the ENTIRE FlowPath tree and fread() one CSV per patient. Before this,
# a full site build did that five times over:
#   build_master()             load_ihc_data()       (result cached as the ihc_data intermediate)
#   report 02                  load_ihc_data()       (uncached)
#   report 04 x2               load_ihc_celltypes()  (uncached — once for ct_all, once for imm_all)
#   report 05 x2               load_ihc_celltypes()  (uncached — same again)
# Reports 04 and 05 each walked the tree TWICE inside a single knit, because ct_all and
# imm_all are separate calls in separate chunks.
#
# Three layers, cheapest first. None changes a call site: the signatures are unchanged and
# every caller keeps working untouched.
#
#   1. SESSION MEMO. A second call in the same R session returns the first call's value.
#      This alone removes the within-report duplication.
#   2. DISK CACHE, staleness-checked. The parsed table is written to output/clean_data/ and
#      re-read on later knits. It is keyed on a STAMP of the source tree — file count, newest
#      mtime, total size — so new or re-exported CSVs invalidate it automatically. A cache
#      that cannot go stale is the failure attend_io.R already has (a fresh .parquet skipped
#      for a months-old .csv), and it is not worth repeating for a slow read.
#   3. COLUMN SELECTION. fread() read every column of a per-cell table; only the handful in
#      FLOWPATH_COLS are ever touched. The header is read first (nrows = 0, effectively free)
#      so only columns that actually exist are requested — the guards downstream still see a
#      genuinely absent marker column as absent.
#
# Escape hatch: pass refresh = TRUE to either loader to force a re-read and rewrite.
# =============================================================================

# Every column either loader touches. phenotype_clean is DERIVED, not read.
FLOWPATH_COLS <- c("phenotype", "Out_of_annotation", "cell_type",
                   "CD45_sign", "CD3_sign", "PD1_sign", "PDL1_sign", "ARID1A_sign")

# Read one FlowPath CSV, requesting only the columns that exist in it. The header probe is
# one line of I/O; the saving is every unused per-cell column across every file.
read_flowpath_csv <- function(csv_path, cols = FLOWPATH_COLS) {
  have <- tryCatch(names(fread(csv_path, nrows = 0L)), error = function(e) NULL)
  keep <- if (is.null(have)) NULL else intersect(cols, have)
  if (length(keep)) fread(csv_path, select = keep) else fread(csv_path)
}

# Cheap fingerprint of the source tree: count, newest mtime, total bytes. Stat calls only —
# no file is opened. NA when the tree is missing, which disables caching rather than
# pretending an empty tree matches.
.flowpath_stamp <- function(dir = FLOWPATH_DIR) {
  if (!dir_exists(dir)) return(NA_character_)
  f <- tryCatch(dir_ls(dir, recurse = TRUE, type = "file", glob = "*.csv"),
                error = function(e) character())
  if (!length(f)) return(NA_character_)
  info <- file.info(f)
  paste0(length(f), "|", as.numeric(max(info$mtime, na.rm = TRUE)), "|",
         sum(info$size, na.rm = TRUE))
}

.ihc_memo <- new.env(parent = emptyenv())

# Session memo -> disk cache -> build. Falls straight through to `builder` when attend_io.R
# has not been sourced (exists() at CALL time, the same rule as maf_standard_cols) or when
# the tree is unreadable, so a bootstrap environment still works.
.flowpath_cached <- function(name, builder, dir = FLOWPATH_DIR, refresh = FALSE) {
  stamp <- .flowpath_stamp(dir)
  memo  <- .ihc_memo[[name]]
  if (!refresh && !is.null(memo) && identical(memo$stamp, stamp)) return(memo$data)

  io_ok <- exists("load_checkpoint", mode = "function") &&
           exists("save_checkpoint", mode = "function")
  if (!refresh && io_ok && !is.na(stamp)) {
    hit <- tryCatch(load_checkpoint(name), error = function(e) NULL)
    if (!is.null(hit) && identical(hit$stamp, stamp) && !is.null(hit$data)) {
      message("load_phenotypes: reusing cached ", name, " (source tree unchanged)")
      .ihc_memo[[name]] <- hit
      return(hit$data)
    }
    if (!is.null(hit)) message("load_phenotypes: ", name, " cache is stale — re-reading.")
  }

  out <- builder()
  obj <- list(stamp = stamp, data = out)
  .ihc_memo[[name]] <- obj
  if (io_ok && !is.na(stamp))
    tryCatch(save_checkpoint(obj, name), error = function(e)
      message("load_phenotypes: could not write ", name, " cache: ", conditionMessage(e)))
  out
}

process_patient <- function(dir) {
  csv_path <- dir_ls(dir, glob = "*.csv")
  if (length(csv_path) != 1)
    stop(sprintf("found %d csv files (expected 1)", length(csv_path)))  # -> skipped by possibly()

  patient_id <- path_file(path_dir(csv_path))
  message(sprintf("LOADING PATIENT: %s", patient_id))

  ihc <- read_flowpath_csv(csv_path) |>
    as_tibble() |>
    mutate(phenotype_clean = str_extract(phenotype, "(?<=\\().*?(?=\\))"))

  inside <- ihc |> filter(Out_of_annotation == FALSE)   # cells within the annotated region

  # PD-1 / PD-L1 staining is only present in some batches, and they are TWO markers with
  # two columns — so they get two guards.
  #
  # ⚠️ There used to be one, `any(str_detect(colnames(inside), "PD1"))`, used to gate BOTH.
  # "PDL1_sign" does not contain the substring "PD1" (P-D-L-1 vs P-D-1), so the guard was
  # driven entirely by the PD-1 column: a batch stained for PD-L1 but not PD-1 had every
  # PD-L1 count set to NA with the data sitting right there. Anchored to the exact column
  # names rather than a substring, so neither can shadow the other again.
  has_pd1  <- "PD1_sign"  %in% colnames(inside)
  has_pdl1 <- "PDL1_sign" %in% colnames(inside)
  # CD3 is part of the same multiplex panel but, like PD-1, absent in some batches.
  has_cd3 <- "CD3_sign" %in% colnames(inside)

  tibble(
    image_id              = patient_id,
    n_cells               = nrow(ihc),
    n_inside              = nrow(inside),
    n_tumor_total         = sum(str_detect(ihc$phenotype, "Tumor")),
    n_tumor_inside        = sum(str_detect(inside$phenotype, "Tumor")),
    # ⚠️ na.rm ON EVERY MARKER COUNT. Without it a single cell whose sign is NA — or whose
    # phenotype carries no parenthesised label, which makes phenotype_clean NA and
    # str_detect() return NA rather than FALSE — turned the WHOLE patient's count into NA.
    # possibly() only catches errors, so that patient then vanished from the panels built on
    # these counts with nothing said. The guards above already distinguish "marker not
    # stained in this batch" (the whole column absent -> NA here) from "this cell was not
    # scored" (dropped from the numerator), which is the distinction that matters.
    n_cd45_inside         = sum(inside$CD45_sign == "+", na.rm = TRUE),
    # CD3+ leukocytes and CD3+CD45+ double-positive T cells inside the annotation.
    n_cd3_inside          = if (has_cd3) sum(inside$CD3_sign == "+", na.rm = TRUE) else NA,
    n_cd3cd45_inside      = if (has_cd3) sum(inside$CD3_sign == "+" & inside$CD45_sign == "+",
                                             na.rm = TRUE) else NA,
    n_tumor_aridia_inside = sum(inside$ARID1A_sign == "+" &
                                  str_detect(inside$phenotype_clean, "Tumor"), na.rm = TRUE),
    # ALL PD-L1-positive cells inside the annotation, whatever their phenotype — the
    # counterpart of n_cd45_inside, and the numerator the tissue-content panel needs. The
    # two counts below it are the tumour and leukocyte SUBSETS of this one; they are kept
    # because the concordance section scores them against the pathologist separately.
    n_pdl1_inside         = if (has_pdl1) sum(inside$PDL1_sign == "+", na.rm = TRUE) else NA,
    n_tumor_pdl1_inside   = if (has_pdl1) sum(inside$PDL1_sign == "+" &
                                               str_detect(inside$phenotype_clean, "Tumor"),
                                             na.rm = TRUE) else NA,
    n_cd45_pd1_inside     = if (has_pd1)  sum(inside$PD1_sign  == "+" & inside$CD45_sign == "+",
                                             na.rm = TRUE) else NA,
    n_cd45_pdl1_inside    = if (has_pdl1) sum(inside$PDL1_sign == "+" & inside$CD45_sign == "+",
                                             na.rm = TRUE) else NA
  )
}

load_ihc_data <- function(flowpath_dir = FLOWPATH_DIR, refresh = FALSE) {
  .flowpath_cached("ihc_data_raw", dir = flowpath_dir, refresh = refresh, builder = function() {
    safe_process <- possibly(process_patient, otherwise = NULL, quiet = FALSE)
    dir_ls(flowpath_dir) |>
      map(safe_process) |>
      list_rbind()
  })
}

# --- Per-cell-type composition (inside the annotation) ----------------------
# For report 05: counts every CELL TYPE inside the annotation, per image, plus the
# three denominators (all cells / tumour cells / CD45+ cells inside) so downstream
# can normalise. "Cell type" = phenotype_clean (the parenthetical label in
# `phenotype`); switch to `phenotype` for the full label. Returns long rows:
# image_id, cell_type, n_cell, n_inside, n_tumor_inside, n_cd45_inside.
process_patient_celltypes <- function(dir) {
  csv_path <- dir_ls(dir, glob = "*.csv")
  if (length(csv_path) != 1)
    stop(sprintf("found %d csv files (expected 1)", length(csv_path)))
  image_id <- path_file(path_dir(csv_path))

  inside <- read_flowpath_csv(csv_path) |>
    as_tibble() |>
    mutate(phenotype_clean = str_extract(phenotype, "(?<=\\().*?(?=\\))")) |>
    filter(Out_of_annotation == FALSE)

  # Exact column names, not substrings — see the guard note in process_patient(): "PDL1_sign"
  # does not contain "PD1", and a substring guard let one marker decide the other's fate.
  has_cd3  <- "CD3_sign"  %in% colnames(inside)
  has_pdl1 <- "PDL1_sign" %in% colnames(inside)

  inside |>
    count(cell_type = phenotype_clean, name = "n_cell") |>
    filter(!is.na(cell_type)) |>
    # na.rm throughout, for the same reason as process_patient(): one unscored cell must not
    # void a whole patient's constant and drop them out of every panel built on it.
    mutate(image_id         = image_id,
           n_inside         = nrow(inside),
           n_tumor_inside   = sum(str_detect(inside$phenotype, "Tumor"), na.rm = TRUE),
           n_cd45_inside    = sum(inside$CD45_sign == "+", na.rm = TRUE),
           # ALL PD-L1+ cells inside — the third column of the tissue-content panel, and the
           # numerator for the PD-L1+ series in ihc_immune_metrics(). NA when this batch was
           # not stained for PD-L1, which is a smaller cohort than the other two columns.
           n_pdl1_inside    = if (has_pdl1) sum(inside$PDL1_sign == "+", na.rm = TRUE) else NA,
           # CD3+CD45+ double-positive T cells inside — denominator for the immune metrics.
           n_cd3cd45_inside = if (has_cd3) sum(inside$CD3_sign == "+" & inside$CD45_sign == "+",
                                               na.rm = TRUE) else NA)
}

load_ihc_celltypes <- function(flowpath_dir = FLOWPATH_DIR, refresh = FALSE) {
  .flowpath_cached("ihc_celltypes", dir = flowpath_dir, refresh = refresh, builder = function() {
    safe_ct <- possibly(process_patient_celltypes, otherwise = NULL, quiet = FALSE)
    dir_ls(flowpath_dir) |>
      map(safe_ct) |>
      list_rbind()
  })
}
