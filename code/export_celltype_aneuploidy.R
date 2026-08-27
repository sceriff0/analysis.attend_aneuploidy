#!/usr/bin/env Rscript
# =============================================================================
# export_celltype_aneuploidy.R  —  the arcsin(sqrt) composition table, as a CSV
#
# One row per patient: pid, aneuploidy_class, then one column per cell type
# holding arcsin(sqrt(fraction of all cells inside the annotation)) — the exact
# quantity report 05's `box-inside-arcsin` panels plot, and report 06's.
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
# Writes output/clean_data/celltype_aneuploidy.csv
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
# annotated region only. frac_inside is the count over every cell in the
# annotation — the "all cells inside" normalisation.
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
wide <- ct |>
  transmute(pid, aneuploidy_class, cell_type,
            asin_frac_inside = arcsin_sqrt(frac_inside)) |>
  pivot_wider(id_cols = c(pid, aneuploidy_class),
              names_from = cell_type, values_from = asin_frac_inside) |>
  arrange(pid)

# Alphabetical cell-type columns, so the header order is stable between runs
# rather than following whichever patient FlowPath happened to read first.
wide <- wide[, c("pid", "aneuploidy_class",
                 sort(setdiff(names(wide), c("pid", "aneuploidy_class"))))]

f_out <- file.path(out_dir, "celltype_aneuploidy.csv")
write_csv(wide, f_out)

cat("\npatients :", nrow(wide), if (mmrd_only) "(MMR-deficient only)" else "", "\n")
cat("cell types:", ncol(wide) - 2L, "\n")
print(table(aneuploidy = wide$aneuploidy_class, useNA = "ifany"))
n_na <- sum(is.na(wide[, -(1:2), drop = FALSE]))
if (n_na > 0)
  cat("\nNOTE:", n_na, "empty patient x cell-type cells — see the comment in this script\n",
      "      before treating them as zeros.\n")
cat("\nwrote:", f_out, "\n")
