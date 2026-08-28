#!/usr/bin/env Rscript
# =============================================================================
# export_celltype_aneuploidy.R  —  the arcsin(sqrt) composition table, as a CSV
#
# One row per patient: pid, aneuploidy_class, then one column per cell type
# holding arcsin(sqrt(n_cell / all cells inside the annotation)) — the exact
# quantity report 05's `box-inside-arcsin` panels plot, and report 06's.
#
# Two extra columns carry tumour and leukocyte content on the SAME denominator:
#   Tumor_cells    n_tumor_inside / n_inside
#   CD45pos_cells  n_cd45_inside  / n_inside
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

args      <- commandArgs(trailingOnly = TRUE)
mmrd_only <- "--mmrd" %in% args
out_dir   <- if ("--out" %in% args) args[which(args == "--out") + 1L] else
               here("output", "clean_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master <- get_master() |> add_molecular_classes()

# Pools each patient's images (counts and denominators summed), inside the
# annotated region only. frac_inside is n_cell / n_inside — the "all cells
# inside" normalisation, the only one exported here.
ct <- ihc_celltype_metrics(load_ihc_celltypes(), load_imaging_data()) |>
  inner_join(master |> select(pid, aneuploidy_class), by = "pid") |>
  filter(!is.na(aneuploidy_class))

if (mmrd_only) {
  # Same is_mmrd() call reports 05/07/09 use, so the restriction cannot drift.
  mmrd_pids <- add_scna_group(master) |> filter(mmr_group == "MMRd") |> pull(pid)
  ct <- ct |> filter(pid %in% mmrd_pids)
}

# --- the per-cell-type columns ----------------------------------------------
# Headers are the cell-type labels verbatim (phenotype_clean, the parenthetical
# label in FlowPath's `phenotype`), so a header matches its facet in the figure.
# A patient with no row for a cell type comes through empty rather than 0: the
# multiplex panel differs between batches, so an absent label can mean "no cells
# of this type" OR "this marker was not stained for this patient", and only the
# FlowPath export can tell those apart.
wide <- ct |>
  transmute(pid, aneuploidy_class, cell_type,
            y = arcsin_sqrt(frac_inside)) |>
  pivot_wider(id_cols = c(pid, aneuploidy_class),
              names_from = cell_type, values_from = y)

# --- tumour and leukocyte content, same denominator -------------------------
# These are NOT cell-type columns and do not partition with them: a tumour cell
# is also counted in whatever phenotype column it belongs to, and the two are
# derived differently upstream — n_tumor_inside matches "Tumor" in `phenotype`,
# n_cd45_inside is CD45_sign == "+" (load_phenotypes.R). They are read off the
# per-patient constants ihc_celltype_metrics() already carries, so they use the
# same pooled n_inside as every column beside them.
extra <- ct |>
  distinct(pid, n_inside, n_tumor_inside, n_cd45_inside) |>
  transmute(pid,
            Tumor_cells   = arcsin_sqrt(n_tumor_inside / n_inside),
            CD45pos_cells = arcsin_sqrt(n_cd45_inside  / n_inside))

clash <- intersect(c("Tumor_cells", "CD45pos_cells"), names(wide))
if (length(clash))
  stop("a FlowPath cell-type label collides with a computed column: ",
       paste(clash, collapse = ", "), " — rename the computed column")

# Tumour and CD45+ lead, then the cell types alphabetically, so the header order
# is stable between runs rather than following whichever patient FlowPath read first.
wide <- wide |>
  left_join(extra, by = "pid") |>
  arrange(pid)
wide <- wide[, c("pid", "aneuploidy_class", "Tumor_cells", "CD45pos_cells",
                 sort(setdiff(names(wide), c("pid", "aneuploidy_class",
                                             "Tumor_cells", "CD45pos_cells"))))]

f_out <- file.path(out_dir, "celltype_aneuploidy.csv")
write_csv(wide, f_out)

cat("\npatients :", nrow(wide), if (mmrd_only) "(MMR-deficient only)" else "", "\n")
cat("columns  :", ncol(wide) - 2L, "(2 computed + ", ncol(wide) - 4L, " cell types)\n")
print(table(aneuploidy = wide$aneuploidy_class, useNA = "ifany"))
n_na <- sum(is.na(as.matrix(wide[, -(1:2), drop = FALSE])))
if (n_na > 0)
  cat("\n", n_na, "empty cells — these are NOT zeros; see the comment above `wide`.\n")
cat("\nwrote:", f_out, "\n")
