#!/usr/bin/env Rscript
# =============================================================================
# export_celltype_aneuploidy.R  —  one CSV: cell-type composition x aneuploidy class
#
# Emits, one row per patient, the arcsine-square-root-transformed cell-type
# fractions that reports 05 (Part 2) and 06 plot, next to that patient's
# aneuploidy class. Nothing is re-derived here: the composition comes from
# ihc_celltype_metrics() and the class from get_master() |> add_molecular_classes(),
# so this export cannot drift from the knitted figures.
#
# Run on the cluster (needs the FlowPath tree + the master's HPC inputs):
#   Rscript code/export_celltype_aneuploidy.R
#   Rscript code/export_celltype_aneuploidy.R --all           # every cell type
#   Rscript code/export_celltype_aneuploidy.R --out /path/dir
#
# Writes to output/clean_data/:
#   celltype_aneuploidy_wide.csv  <- one row per patient  (the one you asked for)
#   celltype_aneuploidy_long.csv  <- one row per patient x cell type (audit trail)
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

source(here("code", "build_master.R"))     # get_master(), add_molecular_classes()
source(here("code", "attend_ihc.R"))       # ihc_celltype_metrics(), arcsin_sqrt()
source(here("code", "load_phenotypes.R"))  # load_ihc_celltypes()
# load_imaging_data() arrives via build_master.R -> load_clinical.R

# --- the cell types wanted, as regexes ---------------------------------------
# "Cell type" is phenotype_clean, the parenthetical label in FlowPath's `phenotype`
# column, and the labels are Italian. Matching is case-insensitive and by pattern
# rather than by literal string because the exact spelling lives in the CSVs on the
# cluster, not in this repo — a literal that misses would drop the column silently.
WANTED <- c(
  "M1"           = "(^|[^A-Za-z0-9])M1([^A-Za-z0-9]|$)",
  "Macrofagi"    = "macrofag",
  "T_citotossici" = "citotoss|cytotox",
  "T_helper"     = "helper"
)

args     <- commandArgs(trailingOnly = TRUE)
keep_all <- "--all" %in% args
out_dir  <- if ("--out" %in% args) args[which(args == "--out") + 1L] else
              here("output", "clean_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. aneuploidy class per patient -----------------------------------------
master <- get_master() |> add_molecular_classes()
# Guarded the way the reports guard it: a missing column yields all-NA rather
# than killing the run, so this stays usable before every attend_cols entry is set.
numcol <- function(df, col) if (col %in% names(df))
  suppressWarnings(as.numeric(df[[col]])) else NA_real_
aneu <- tibble(pid              = master$pid,
               aneuploidy_class = master$aneuploidy_class,
               aneuploidy_score = numcol(master, attend_cols$aneuploidy))

# --- 2. cell-type composition per patient ------------------------------------
# Pools each patient's images (counts and denominators summed), inside the
# annotated region only. Three denominators, because each asks a different
# question — see report 05 Part 2.
ct <- ihc_celltype_metrics(load_ihc_celltypes(), load_imaging_data())

cat("\ncell types observed in the FlowPath export:\n")
print(sort(unique(ct$cell_type)))

# --- 3. restrict to the requested cell types ---------------------------------
if (!keep_all) {
  hits <- imap_dfr(WANTED, function(pat, label) {
    m <- unique(ct$cell_type[str_detect(ct$cell_type, regex(pat, ignore_case = TRUE))])
    if (length(m) == 0) {
      warning("no cell type matched '", label, "' (pattern: ", pat, ")", call. = FALSE)
      return(tibble())
    }
    tibble(cell_type = m, label = label)
  })
  cat("\nrequested label -> matched cell type(s):\n")
  if (nrow(hits)) print(as.data.frame(hits)) else cat("  (nothing matched)\n")
  # A cell type matching two labels would be counted twice downstream; flag it.
  dup <- hits$cell_type[duplicated(hits$cell_type)]
  if (length(dup)) warning("cell type(s) matched by more than one label: ",
                           paste(unique(dup), collapse = ", "), call. = FALSE)
  ct <- ct |> inner_join(hits, by = "cell_type")
} else {
  ct <- ct |> mutate(label = cell_type)
}

# --- 4. long table: raw + arcsine fractions ----------------------------------
long <- ct |>
  inner_join(aneu, by = "pid") |>
  transmute(pid, cell_type, label, aneuploidy_class, aneuploidy_score,
            n_cell, n_inside, n_tumor_inside, n_cd45_inside,
            frac_inside, frac_tumor, frac_cd45,
            asin_frac_inside = arcsin_sqrt(frac_inside),
            asin_frac_tumor  = arcsin_sqrt(frac_tumor),
            asin_frac_cd45   = arcsin_sqrt(frac_cd45)) |>
  arrange(pid, label)

# --- 5. wide table: one row per patient --------------------------------------
# Column names are <label>__<metric>, e.g. M1__asin_frac_inside. make.names() is
# NOT used on the keyed columns (pid, aneuploidy_*) — only the cell-type labels
# are sanitised, and the master's prefixed names never enter this table.
wide <- long |>
  select(pid, aneuploidy_class, aneuploidy_score, label,
         frac_inside, frac_tumor, frac_cd45,
         asin_frac_inside, asin_frac_tumor, asin_frac_cd45) |>
  mutate(label = make.names(label)) |>
  pivot_wider(id_cols = c(pid, aneuploidy_class, aneuploidy_score),
              names_from = label,
              values_from = c(frac_inside, frac_tumor, frac_cd45,
                              asin_frac_inside, asin_frac_tumor, asin_frac_cd45),
              names_glue = "{label}__{.value}") |>
  arrange(pid)

f_long <- file.path(out_dir, "celltype_aneuploidy_long.csv")
f_wide <- file.path(out_dir, "celltype_aneuploidy_wide.csv")
write_csv(long, f_long)
write_csv(wide, f_wide)

cat("\npatients with both a cell-type measurement and an aneuploidy class:",
    n_distinct(long$pid), "\n")
print(table(aneuploidy = long |> distinct(pid, aneuploidy_class) |>
              pull(aneuploidy_class), useNA = "ifany"))
cat("\nwrote:\n  ", f_wide, "\n  ", f_long, "\n")
