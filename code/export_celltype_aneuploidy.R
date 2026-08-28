#!/usr/bin/env Rscript
# =============================================================================
# export_celltype_aneuploidy.R  —  the arcsin(sqrt) composition table, as a CSV
#
# One row per patient: pid, aneuploidy_class, then one column per cell type
# holding arcsin(sqrt(fraction)) — the exact quantity report 05's `box-*-arcsin`
# panels plot, and report 06's. One file per denominator, because the denominator
# changes the question and a single table mixing them buries that:
#
#   _inside.csv  n_cell / all cells inside     what fraction of the tissue is this
#   _tumor.csv   n_cell / tumour cells inside  this cell type per unit of tumour
#   _cd45.csv    n_cell / CD45+ cells inside   composition OF the immune compartment
#
# Nothing is re-derived: the composition comes from ihc_celltype_metrics() and
# the class from get_master() |> add_molecular_classes(), so this CSV cannot
# drift from the knitted figures.
#
# Run on the cluster (needs the FlowPath tree + the master's HPC inputs):
#   Rscript code/export_celltype_aneuploidy.R
#   Rscript code/export_celltype_aneuploidy.R --mmrd        # MMRd only, as in report 05
#   Rscript code/export_celltype_aneuploidy.R --out /path/dir
#
# Writes output/clean_data/celltype_aneuploidy_{inside,tumor,cd45}.csv
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

source(here("code", "build_master.R"))     # get_master(), add_molecular_classes(), add_scna_group()
source(here("code", "attend_ihc.R"))       # ihc_celltype_metrics(), arcsin_sqrt()
source(here("code", "load_phenotypes.R"))  # load_ihc_celltypes()
# load_imaging_data() arrives via build_master.R -> load_clinical.R

args     <- commandArgs(trailingOnly = TRUE)
mmrd_only <- "--mmrd" %in% args
out_dir  <- if ("--out" %in% args) args[which(args == "--out") + 1L] else
              here("output", "clean_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master <- get_master() |> add_molecular_classes()

# Pools each patient's images (counts and denominators summed), inside the
# annotated region only. Carries all three fractions — frac_inside / frac_tumor /
# frac_cd45 — so the three exports below share one count and one pooling step.
ct <- ihc_celltype_metrics(load_ihc_celltypes(), load_imaging_data()) |>
  inner_join(master |> select(pid, aneuploidy_class), by = "pid") |>
  filter(!is.na(aneuploidy_class))

if (mmrd_only) {
  # Same is_mmrd() call reports 05/07/09 use, so the restriction cannot drift.
  mmrd_pids <- add_scna_group(master) |> filter(mmr_group == "MMRd") |> pull(pid)
  ct <- ct |> filter(pid %in% mmrd_pids)
}

# Columns are the cell-type labels verbatim (phenotype_clean, the parenthetical
# label in FlowPath's `phenotype`), so a header matches its facet in the figure.
# A patient with no row for a cell type comes through NA rather than 0: the
# multiplex panel differs between batches, so an absent label can mean "no cells
# of this type" OR "this marker was not stained for this patient", and only the
# FlowPath export can tell those apart.
#
# arcsin_sqrt() clamps to [0, 1] before asin(sqrt()). That is not cosmetic for the
# tumour and CD45+ denominators: a cell type more abundant than the denominator
# gives a fraction above 1, which asin() would return NaN for. Clamped, it lands
# on pi/2 — a visible ceiling in the data rather than a silent hole. The run
# reports how many cells hit it.
write_one <- function(frac_col, suffix) {
  wide <- ct |>
    transmute(pid, aneuploidy_class, cell_type,
              y = arcsin_sqrt(.data[[frac_col]])) |>
    pivot_wider(id_cols = c(pid, aneuploidy_class),
                names_from = cell_type, values_from = y) |>
    arrange(pid)

  # Alphabetical cell-type columns, so the header order is stable between runs
  # rather than following whichever patient FlowPath happened to read first.
  wide <- wide[, c("pid", "aneuploidy_class",
                   sort(setdiff(names(wide), c("pid", "aneuploidy_class"))))]

  f_out <- file.path(out_dir, paste0("celltype_aneuploidy_", suffix, ".csv"))
  write_csv(wide, f_out)

  vals    <- as.matrix(wide[, -(1:2), drop = FALSE])
  n_na    <- sum(is.na(vals))
  n_ceil  <- sum(!is.na(vals) & vals >= pi / 2 - 1e-9)
  cat(sprintf("%-8s %3d patients x %2d cell types  |  %d empty  |  %d clamped at pi/2  ->  %s\n",
              suffix, nrow(wide), ncol(wide) - 2L, n_na, n_ceil, basename(f_out)))
  invisible(wide)
}

cat("\npatients:", n_distinct(ct$pid), if (mmrd_only) "(MMR-deficient only)" else "", "\n")
print(table(aneuploidy = ct |> distinct(pid, aneuploidy_class) |> pull(aneuploidy_class),
            useNA = "ifany"))
cat("\n")

# (fraction column on the ct table, file suffix) — the three denominators of
# report 05 Part 2, in the order that report shows them.
invisible(purrr::imap(c(inside = "frac_inside", tumor = "frac_tumor", cd45 = "frac_cd45"),
                      function(col, suffix) write_one(col, suffix)))

cat("\nEmpty cells are NOT zeros — see the comment above write_one() before filling them.\n")
cat("wrote to:", out_dir, "\n")
