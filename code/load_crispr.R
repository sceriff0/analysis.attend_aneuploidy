# =============================================================================
# load_crispr.R — the Wu et al. (Immunity 2025) resister / sensitizer catalogue
#
# Exports load_crispr_screens(). Reads the tracked TSV that
# code/fetch_crispr_screens.R derived from the paper's supplementary workbook
# (Table S3, 1-s2.0-S1074761325004315-mmc4.xlsx) and applies the measurement cutoff as a
# DECLARED argument.
#
# Returns an empty frame with the right columns when the file is absent, so every report
# knits before the data lands — the same contract as the load_*.R loaders (CLAUDE.md,
# "Missing data-loading layer").
#
# Config: attend_crispr (code/attend_classes.R). Design:
# specs/2026-09-24-wu2025-crispr-functional-filter.md.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(readr)
})

#' The CRISPR T-cell-killing catalogue, filtered to the configured evidence level.
#'
#' @param min_measurements keep genes with at least this many concordant independent
#'   measurements. Wu et al.'s own cutoff is 2, which reproduces their reported 519
#'   resisters / 877 sensitizers exactly. Passed explicitly so a sensitivity check is one
#'   argument, not an edit to a tracked data file.
#' @param drop_both_list drop gene symbols that qualify as BOTH resister and sensitizer.
#'   They are still counted and returned via attr(, "both_list") so the report can print
#'   how many were dropped — an unexplained absence is the thing note_skip() exists for.
#' @param cfg attend_crispr.
#'
#' @return tibble(gene, role, n_measurements, wu_overlap) where `role` is
#'   "resister"/"sensitizer" and `wu_overlap` is Wu et al.'s OWN result on their melanoma
#'   cohort. Attributes: `both_list` (the dropped symbols), `n_unfiltered`, `min_measurements`.
load_crispr_screens <- function(min_measurements = attend_crispr$min_measurements,
                                drop_both_list   = attend_crispr$drop_both_list,
                                cfg              = attend_crispr) {
  empty <- tibble(gene = character(0), role = character(0),
                  n_measurements = integer(0), wu_overlap = character(0))

  path <- here("data", cfg$file)
  if (!file.exists(path)) {
    message("load_crispr_screens(): no ", path, " — run Rscript code/fetch_crispr_screens.R")
    return(empty)
  }

  raw <- tryCatch(
    read_tsv(path, col_types = cols(gene = col_character(), role = col_character(),
                                    n_measurements = col_integer(),
                                    wu2025_dp_overlap = col_character()),
             progress = FALSE),
    error = function(e) { message("load_crispr_screens(): ", conditionMessage(e)); NULL })
  if (is.null(raw) || !nrow(raw)) return(empty)

  # ⚠️ THE SOURCE COLUMN IS TWO VARIABLES UNDER ONE NAME. In the workbook, column M is
  # "Overlap status with recurrent, DP-specific, DELETED genes" (the resister block) and
  # column T is "...AMPLIFIED genes" (the sensitizer block). The converter writes both into
  # wu2025_dp_overlap, so the name alone does not say which direction a row refers to —
  # `role` does. Renamed here to wu_overlap and documented, rather than left looking like
  # one variable. This is the shape of the .peak_altered() bug: a direction column that was
  # really a group label under the wrong name.
  out <- raw |>
    rename(wu_overlap = "wu2025_dp_overlap") |>
    filter(!is.na(.data$n_measurements))

  n_unfiltered <- nrow(out)
  out <- filter(out, .data$n_measurements >= as.integer(min_measurements))

  # A symbol reported as a resister by some screens and a sensitizer by others cannot
  # support a directional claim, and Wu et al. give no rule for them. Same treatment as
  # gistic_feature_direction()'s amp-and-del genes: excluded, never silently resolved by
  # whichever row a join happened to keep.
  both <- intersect(out$gene[out$role == "resister"], out$gene[out$role == "sensitizer"])
  if (drop_both_list && length(both)) out <- filter(out, !.data$gene %in% both)

  attr(out, "both_list")        <- both
  attr(out, "n_unfiltered")     <- n_unfiltered
  attr(out, "min_measurements") <- as.integer(min_measurements)
  out
}

#' Counts of the catalogue as a printable partition — the companion counter.
#'
#' Prints beside every figure that uses the catalogue, so the reader sees how many genes
#' the evidence cutoff removed and how many were ambiguous, rather than one filtered number
#' with no denominator. Same role as mutation_counts() / response_counts() / promise_counts().
crispr_counts <- function(crispr = load_crispr_screens()) {
  both <- attr(crispr, "both_list")
  data.frame(
    quantity = c(paste0("catalogue rows, unfiltered"),
                 paste0("kept at >= ", attr(crispr, "min_measurements"), " measurements"),
                 "  of which resisters (lost -> immune escape)",
                 "  of which sensitizers (gained -> immune escape)",
                 "dropped: symbol in BOTH lists, no direction"),
    n = c(as.integer(attr(crispr, "n_unfiltered")),
          nrow(crispr),
          sum(crispr$role == "resister"),
          sum(crispr$role == "sensitizer"),
          length(both)),
    stringsAsFactors = FALSE
  )
}
