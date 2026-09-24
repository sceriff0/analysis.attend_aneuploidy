#!/usr/bin/env Rscript
# =============================================================================
# fetch_crispr_screens.R — Wu et al. (Immunity 2025) Table S3 -> a tracked TSV
#
# Converts the supplementary workbook (1-s2.0-S1074761325004315-mmc4.xlsx) into
# data/crispr_t_cell_screens_wu2025.tsv. Run ONCE, per the same pattern as
# code/fetch_tcga_ucec_2013.R: the .xlsx is raw supplementary material and stays out of
# git; the small derived TSV is a reference table and IS tracked, exactly as
# data/tcga_gistic_peaks_2013.tsv is.
#
# WHAT THE TABLE IS. Wu et al. pooled 23 published CRISPR-Cas9 screens of tumour-cell
# killing by T cells and classified every hit as a RESISTER gene (knockout promoted
# resistance) or a SENSITIZER gene (knockout promoted sensitivity). Their recurrent-CNV
# result is the intersection of recurrently DELETED genes with resisters, and recurrently
# AMPLIFIED genes with sensitizers — the functional filter that turns a several-thousand
# gene recurrence list into a finding.
#
# ⚠️ THE WORKBOOK SHIPS UNFILTERED, AND THE UNFILTERED SET IS 5x TOO BIG. The sheet
# carries 2,828 resister and 2,961 sensitizer rows, with `Number of measurements` as low
# as 1. The paper's 519 / 877 appear only after its stated cutoff — "concordant evidence
# in >= two independent measurements" (STAR Methods, "Cataloging resister and sensitizer
# genes"). Verified here: >= 2 reproduces 519 and 877 exactly. The cutoff is therefore NOT
# applied in this converter — the full table is written with n_measurements intact and the
# filter is a DECLARED argument at load time (attend_crispr$min_measurements), so the
# choice is visible in config rather than baked into a data file nobody re-reads.
#
# ⚠️ 105 GENES ARE IN BOTH LISTS at that cutoff — reported as resisters by some screens and
# sensitizers by others. The paper allows this (519 + 877 do not deduplicate), but a gene
# that is both cannot support a DIRECTIONAL claim, which is the only kind this analysis
# makes. They are written here and marked downstream, never silently resolved by whichever
# row a join happened to keep.
#
# The `wu2025_dp_overlap` column is Wu et al.'s OWN result on their melanoma cohort
# (YES = the gene was recurrently, DP-specifically deleted/amplified in their patients).
# It is carried for provenance and as a regression anchor — 108 YES on the resister side,
# matching the paper's Figure 1D — and it says NOTHING about ATTEND. Never join on it.
#
# RUN (repo root):  Rscript code/fetch_crispr_screens.R [path/to/mmc4.xlsx]
# =============================================================================
suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(tibble)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
xlsx <- if (length(args) >= 1) args[[1]] else {
  # The publisher's filename is stable but ugly; glob it rather than hardcode the DOI stem.
  hits <- Sys.glob(here("*mmc4.xlsx"))
  if (length(hits) == 0)
    stop("no *mmc4.xlsx in ", here(), " — pass the path as an argument.")
  hits[[1]]
}
if (!file.exists(xlsx)) stop("missing workbook: ", xlsx)

out <- here("data", "crispr_t_cell_screens_wu2025.tsv")

# The sheet is THREE blocks side by side in one grid, not a rectangle: cols A-G are the
# study metadata, J-M the resister genes, Q-T the sensitizer genes. Row 1 holds the block
# titles and row 2 the real header, so the data starts at row 3. Reading with col_names =
# FALSE and skip = 2 keeps the positional layout — naming the columns from row 2 would
# collide, since both gene blocks call their columns "Gene".
raw <- read_excel(xlsx, sheet = 1, col_names = FALSE, skip = 2,
                  .name_repair = "minimal")

# Positional extraction of one gene block. `cols` are 1-based column indices into the grid:
# gene, number of measurements, publication source details, overlap status.
block <- function(dat, cols, role) {
  g <- as.character(dat[[cols[1]]])
  keep <- !is.na(g) & nzchar(trimws(g))
  tibble(
    gene              = trimws(g[keep]),
    role              = role,
    n_measurements    = suppressWarnings(as.integer(as.character(dat[[cols[2]]])[keep])),
    wu2025_dp_overlap = trimws(as.character(dat[[cols[3]]])[keep])
  )
}

# J,K,M = 10,11,13   Q,R,T = 17,18,20.  Column L / S (the free-text publication details)
# is deliberately DROPPED: up to 1.2 kB of prose per row, with full-width colons as its
# internal separator, and nothing downstream reads it. The PMIDs in block A-G are the
# citable provenance; the per-gene prose is not.
screens <- bind_rows(
  block(raw, c(10, 11, 13), "resister"),
  block(raw, c(17, 18, 20), "sensitizer")
) |>
  filter(!is.na(n_measurements)) |>
  arrange(role, desc(n_measurements), gene)

write_tsv(screens, out)

# ---- verification against the paper, printed every run -----------------------
# These four numbers are stated in the paper (519 / 877 in the Results; 108 / 90 in
# Figure 1D-1E). If a re-download or a re-save of the workbook shifts a column, they move,
# and the conversion is wrong in a way the file size would never show.
n_at <- function(role, cut) sum(screens$role == role & screens$n_measurements >= cut)
r2 <- n_at("resister", 2); s2 <- n_at("sensitizer", 2)
both <- length(intersect(
  screens$gene[screens$role == "resister"   & screens$n_measurements >= 2],
  screens$gene[screens$role == "sensitizer" & screens$n_measurements >= 2]))
ov <- function(role) sum(screens$role == role & screens$n_measurements >= 2 &
                           screens$wu2025_dp_overlap == "YES")

cat("wrote ", out, "\n",
    "  rows (unfiltered)          : ", nrow(screens), "\n",
    "  resister    >=2 measurements: ", r2, "   (paper: 519)\n",
    "  sensitizer  >=2 measurements: ", s2, "   (paper: 877)\n",
    "  in BOTH lists at >=2        : ", both, "   <- not directionally interpretable\n",
    "  resister   & DP-deleted     : ", ov("resister"),   "   (paper: 108)\n",
    "  sensitizer & DP-amplified   : ", ov("sensitizer"), "   (paper:  90)\n", sep = "")

if (r2 != 519L || s2 != 877L)
  stop("Table S3 conversion does not reproduce the paper's 519/877 — check the column ",
       "offsets against the workbook before using this file.")

# ⚠️ 519/877 ALONE CANNOT CATCH A DIRECTION SWAP. Those two counts depend only on the
# `Number of measurements` offsets (K and R). The overlap flags live in the ADJACENT columns
# M and T, and they are different variables: M is "Overlap status with recurrent,
# DP-specific, DELETED genes" (resisters), T is "...AMPLIFIED genes" (sensitizers). Swap M
# and T and 519/877 still pass while the direction of every downstream claim inverts —
# exactly the shape of the .peak_altered() sign bug, which the fixtures also failed to catch.
# 108 and 90 are stated in the paper (Figures 1D and 1E), so they pin the orientation too.
if (ov("resister") != 108L || ov("sensitizer") != 90L)
  stop("Table S3 conversion reproduces 519/877 but NOT the paper's 108/90 overlap counts ",
       "(got ", ov("resister"), "/", ov("sensitizer"), ") — columns M and T have most ",
       "likely been swapped, which silently inverts every directional result downstream.")
