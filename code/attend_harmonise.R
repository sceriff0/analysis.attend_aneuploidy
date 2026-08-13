# =============================================================================
# attend_harmonise.R
# Shared harmonisation logic for the ATTEND QC pipeline.
#
# Sourced by BOTH:
#   analysis/01-data-integration.Rmd  (dataset assembly)
#   analysis/01-data-integration.Rmd  (Part 2: linkage diagnostics + content QC)
#
# Keeping these definitions in one place guarantees the two reports use
# identical ID spaces, crosswalks and collapsing rules — they cannot drift
# apart the way two copy-pasted code blocks would.
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))

# --- ID configuration -------------------------------------------------------
# For each dataset: which column holds its native ID, and which ID space that
# ID lives in (pid = patient, barcode = sequencing sample, image = slide).
make_id_cfg <- function() list(
  tmb     = list(col = "ID",            space = "barcode"),
  msi     = list(col = "ID",            space = "barcode"),
  maf     = list(col = "ID",            space = "barcode"),
  hrd     = list(col = "Sample",        space = "barcode"),
  aneu    = list(col = "sequencing_id", space = "barcode"),
  gianlu  = list(col = "ID",            space = "pid",   barcode_col = "TUMOR_BARCODE"),
  imaging = list(col = "ID",            space = "pid",   image_col   = "image_id"),
  ihc     = list(col = "image_id",      space = "image")
)

# The genomic tables keyed by barcode (used by several QC helpers).
barcode_datasets <- c("tmb", "msi", "maf", "hrd", "aneu")

# --- Crosswalk builders -----------------------------------------------------
# barcode -> pid, reverse-engineered from gianlu_clinical_data (which carries
# both TUMOR_BARCODE and the patient ID). NOTE: the canonical map is
# attend_barcodes.csv; 01-data-integration.Rmd validates this derived crosswalk
# against it.
build_barcode_pid <- function(gianlu_clinical_data, cfg) {
  gianlu_clinical_data |>
    transmute(barcode = as.character(.data[[ cfg$gianlu$barcode_col ]]),
              pid     = as.character(.data[[ cfg$gianlu$col ]])) |>
    filter(!is.na(barcode), !is.na(pid)) |>
    distinct()
}

# image -> pid, from imaging_data (carries both image_id and patient ID).
build_image_pid <- function(imaging_data, cfg) {
  imaging_data |>
    transmute(image = as.character(.data[[ cfg$imaging$image_col ]]),
              pid   = as.character(.data[[ cfg$imaging$col ]])) |>
    filter(!is.na(image), !is.na(pid)) |>
    distinct()
}

# --- Normalisation for safe ID recovery -------------------------------------
# IDENTICAL to report 01's near-miss norm() so the recovery and the ledger that
# reports it cannot drift. Upper-cases and trims; nothing else (suffixes are
# already stripped by the loaders before linkage).
norm_id <- function(x) toupper(trimws(as.character(x)))

# Build a normalised crosswalk key -> pid, KEEPING ONLY norm-keys that map to a
# single distinct pid. Ambiguous norm-keys (>1 pid) are dropped = refused, so a
# near-miss can never be mapped to the wrong patient. Returns a named character
# vector (nkey -> pid).
.build_norm_cw <- function(cw, key_col) {
  k   <- norm_id(cw[[key_col]])
  ok  <- !is.na(k) & !is.na(cw$pid)
  k   <- k[ok]; p <- as.character(cw$pid[ok])
  by  <- tapply(p, k, function(v) unique(v))
  uniq <- by[vapply(by, length, integer(1)) == 1]
  stats::setNames(vapply(uniq, `[[`, character(1), 1L), names(uniq))
}

# --- Translator -------------------------------------------------------------
# Returns a closure that maps any dataset's native IDs onto pid. NA = could not
# be linked to a patient. `recover = TRUE` adds a collision-guarded normalised
# second pass for IDs the exact match missed (case/whitespace drift only).
make_pid_vector <- function(cw_barcode_pid, cw_image_pid, recover = FALSE) {
  norm_barcode <- if (recover) .build_norm_cw(cw_barcode_pid, "barcode") else NULL
  norm_image   <- if (recover) .build_norm_cw(cw_image_pid,   "image")   else NULL
  function(df, cfg) {
    ids <- as.character(df[[cfg$col]])
    if (cfg$space == "pid") return(ids)
    # switch() returns NULL for an unrecognised space, exactly as the original did —
    # an unknown/mistyped space must fail, never silently borrow the wrong crosswalk.
    cw   <- switch(cfg$space, barcode = cw_barcode_pid, image = cw_image_pid)
    if (is.null(cw)) return(NULL)
    key  <- if (cfg$space == "barcode") "barcode" else "image"
    out  <- cw$pid[match(ids, cw[[key]])]
    if (recover) {
      miss <- is.na(out)
      if (any(miss)) {
        ncw <- switch(cfg$space, barcode = norm_barcode, image = norm_image)
        out[miss] <- unname(ncw[match(norm_id(ids[miss]), names(ncw))])
      }
    }
    out
  }
}

# Convenience: the patient set reached by each dataset.
build_sets <- function(tables, id_cfg, pid_vector) {
  Map(function(df, cfg) {
    pid <- pid_vector(df, cfg)
    unique(pid[!is.na(pid)])
  }, tables, id_cfg[names(tables)])
}

# --- Collapsing rule --------------------------------------------------------
# When a patient has several rows (multiple samples/images, or many MAF
# variants), numeric columns are averaged and non-numeric take the first
# non-missing value. This is LOSSY by design — 01-data-integration.Rmd Part 2 reports
# exactly how many patients are affected.
num_or_first <- function(x)
  if (is.numeric(x)) mean(x, na.rm = TRUE) else x[which(!is.na(x))][1]

# Translate a dataset to pid, keep mapped rows, collapse to one row per patient,
# and prefix payload columns with the dataset name to avoid collisions.
collapse_pid <- function(df, cfg, prefix, pid_vector) {
  pid  <- pid_vector(df, cfg)
  keep <- !is.na(pid)
  df2  <- df[keep, , drop = FALSE]
  df2$.pid <- pid[keep]
  drop_cols <- c(cfg$col, cfg$barcode_col, cfg$image_col)   # NULLs drop out of c()
  payload   <- setdiff(names(df2), c(drop_cols, ".pid"))
  if (length(payload) == 0) return(tibble(pid = unique(df2$.pid)))
  out <- df2 |>
    group_by(.pid) |>
    summarise(across(all_of(payload), num_or_first), .groups = "drop")
  names(out)[-1] <- paste0(prefix, "__", names(out)[-1])
  rename(out, pid = .pid)
}

# --- Generic numeric content audit ------------------------------------------
# Per numeric column: NA/zero/negative rates, range, distinct count, and a
# robust (median ± 5*MAD) outlier count. Used for value-range sanity checks.
.safe_stat <- function(x, f) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)

numeric_audit <- function(df, dataset) {
  num <- dplyr::select(df, where(is.numeric))
  if (ncol(num) == 0) return(tibble())
  purrr::imap(num, function(x, col) {
    med  <- .safe_stat(x, stats::median)
    madv <- .safe_stat(x, stats::mad)
    tibble(
      dataset       = dataset,
      column        = col,
      n             = length(x),
      pct_na        = round(mean(is.na(x)) * 100, 1),
      n_distinct    = dplyr::n_distinct(x, na.rm = TRUE),
      min           = .safe_stat(x, min),
      median        = med,
      max           = .safe_stat(x, max),
      pct_zero      = round(mean(x == 0, na.rm = TRUE) * 100, 1),
      pct_negative  = round(mean(x < 0,  na.rm = TRUE) * 100, 1),
      n_outlier_mad = if (is.na(madv) || madv == 0) 0L
                      else sum(abs(x - med) > 5 * madv, na.rm = TRUE)
    )
  }) |> bind_rows()
}
