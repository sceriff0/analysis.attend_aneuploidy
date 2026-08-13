# =============================================================================
# attend_io.R  —  format-agnostic on-disk I/O for pipeline intermediates.
#
# Sourced by report 01 (writes) and every report from 11 on (reads).
#
# Prefers Apache Parquet via `arrow` — typed, columnar, compressed (~10x smaller
# than CSV), and immune to the schema drift that bites a CSV re-read across six
# reports (see research/ §6). Falls back to CSV when arrow isn't installed, so
# every report still runs before the dependency is added. Reads prefer .parquet
# and fall back to .csv, so a half-migrated output/ directory keeps working.
#
# Use BASE NAMES (no extension): write_intermediate(df, "attend_master_joined").
# =============================================================================

suppressPackageStartupMessages(library(here))

.attend_clean_dir <- function() here::here("output", "clean_data")

# Write `df` as an intermediate called `name`. Returns the path actually written.
write_intermediate <- function(df, name, dir = .attend_clean_dir()) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  if (requireNamespace("arrow", quietly = TRUE)) {
    path <- file.path(dir, paste0(name, ".parquet"))
    arrow::write_parquet(df, path)
  } else {
    path <- file.path(dir, paste0(name, ".csv"))
    utils::write.csv(df, path, row.names = FALSE)
  }
  invisible(path)
}

# Read intermediate `name`, preferring .parquet then .csv.
read_intermediate <- function(name, dir = .attend_clean_dir()) {
  pq  <- file.path(dir, paste0(name, ".parquet"))
  csv <- file.path(dir, paste0(name, ".csv"))
  if (file.exists(pq) && requireNamespace("arrow", quietly = TRUE)) {
    tibble::as_tibble(arrow::read_parquet(pq))
  } else if (file.exists(csv)) {
    readr::read_csv(csv, show_col_types = FALSE)
  } else {
    stop("No intermediate '", name, "' in ", dir, " (looked for .parquet and .csv). ",
         "Run 01_data_integration.Rmd first.")
  }
}

# R<->R object checkpoints (e.g. handing loaded HPC tables between reports): qs2
# if available (fast saveRDS replacement; the original `qs` is deprecated), else RDS.
save_checkpoint <- function(obj, name, dir = .attend_clean_dir()) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  if (requireNamespace("qs2", quietly = TRUE)) {
    path <- file.path(dir, paste0(name, ".qs2")); qs2::qs_save(obj, path)
  } else {
    path <- file.path(dir, paste0(name, ".rds")); saveRDS(obj, path)
  }
  invisible(path)
}

load_checkpoint <- function(name, dir = .attend_clean_dir()) {
  qs  <- file.path(dir, paste0(name, ".qs2"))
  rds <- file.path(dir, paste0(name, ".rds"))
  if (file.exists(qs) && requireNamespace("qs2", quietly = TRUE)) qs2::qs_read(qs)
  else if (file.exists(rds)) readRDS(rds)
  else stop("No checkpoint '", name, "' in ", dir, ".")
}
