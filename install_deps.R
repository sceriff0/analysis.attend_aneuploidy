#!/usr/bin/env Rscript
# =============================================================================
# install_deps.R — dependency installer for the ATTEND Docker image
# =============================================================================
# Two modes, chosen automatically:
#
#   1. renv.lock present  -> renv::restore()  (reproducible, exact versions)
#   2. renv.lock absent   -> install the curated set below
#
# The curated set is derived directly from what the reports and code/ scripts
# load (library/require/requireNamespace/`pkg::`), including the uncommitted
# changes (notably code/attend_io.R's arrow/qs2 path). Keep this in sync with
# the source until `renv::init(bioconductor = "3.18")` makes renv.lock the
# single source of truth — after that, this list is just a fallback.
# =============================================================================

options(
  repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"),
  Ncpus = max(1L, parallel::detectCores())
)

install_if_missing <- function(pkgs, installer) {
  have <- rownames(installed.packages())
  need <- setdiff(pkgs, have)
  if (length(need)) installer(need)
}

if (file.exists("renv.lock")) {
  message("== renv.lock found -> restoring pinned environment ==")
  if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
  renv::restore(prompt = FALSE)
} else {
  message("== No renv.lock -> installing curated dependency set ==")

  # BiocManager is preconfigured to release 3.18 in the base image.
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

  # --- CRAN ----------------------------------------------------------------
  cran <- c(
    # core + project plumbing
    "tidyverse",      # dplyr / ggplot2 / purrr / readr / tibble ...
    "here", "fs", "renv", "workflowr", "knitr", "rmarkdown",
    # I/O intermediates (code/attend_io.R)
    "arrow", "qs2",
    # report 01 — ID linkage QC
    "UpSetR", "reclin2",
    # report 02 — master dataset QC/EDA
    "naniar", "skimr", "pointblank",
    # reports 03/05/06 — survival
    "survival", "ggsurvfit", "gtsummary", "broom",
    # report 04 — aneuploidy vs classes
    "ggpubr",
    # report 06 — spatial
    "spatstat.geom", "sf", "dbscan"
  )
  install_if_missing(cran, install.packages)

  # --- Bioconductor 3.18 ---------------------------------------------------
  bioc <- c(
    "maftools",          # code/attend_classes.R, load_wes_results.R (real MAFs)
    "SpatialExperiment", # report 05/06 spatial stack
    "imcRtools"          # report 05/06 spatial stack
    # "spicyR"           # noted in CLAUDE.md for future spatial work — add when needed
  )
  install_if_missing(bioc, function(p) BiocManager::install(p, update = FALSE, ask = FALSE))
}

message("== dependency installation complete ==")
