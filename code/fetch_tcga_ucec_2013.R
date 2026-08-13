#!/usr/bin/env Rscript
# =============================================================================
# fetch_tcga_ucec_2013.R — the TCGA UCEC *2013 Nature paper* cohort (report 07)
#
# Sibling of fetch_tcga_ucec.R, which pulls the PanCancer-Atlas (2018) restatement.
# THIS script pulls the ORIGINAL publication study — cBioPortal `ucec_tcga_pub`,
# "Uterine Corpus Endometrial Carcinoma (TCGA, Nature 2013)", PMID 23636398, the exact
# cohort behind papers/nature12113.pdf (373 samples, 248 exomes, reference genome hg19).
#
# Four files land in data/tcga_2013/:
#   data_clinical_sample.txt   CNA_CLUSTER_K4 (the Fig-1a copy-number clusters, n=363),
#                              MSI_STATUS_7_MARKER_CALL, MLH1_SILENCING,
#                              MUTATION_RATE_CLUSTER, FRACTION_GENOME_ALTERED, TMB
#   data_clinical_patient.txt  SUBTYPE = the paper's four integrated groups (n=232):
#                              POLE (Ultra-mutated) 17 / MSI (Hyper-mutated) 65 /
#                              Copy-number low (Endometriod) 90 /
#                              Copy-number high (Serous-like) 60
#   data_cna_hg19.seg          SNP6 segmented copy number, ~52k segments x 365 samples —
#                              the rows (chromosomal location) of Fig. 1a
#   pancan_aneuploidy.tsv      ANEUPLOIDY_SCORE, which the 2013 study does NOT carry.
#
# WHY THE FOURTH FILE. The 2013 paper predates per-sample aneuploidy scoring; the
# canonical arm-level score (0-39) is Taylor et al., Cancer Cell 2018, published on
# cBioPortal only under ucec_tcga_pan_can_atlas_2018. We fetch it there and join on
# PATIENT_ID (TCGA-XX-XXXX), which recovers 225 of the 232 classified patients. It is
# therefore a LATER measure retrofitted onto the 2013 cohort — reports say so, and
# also plots the study-native FRACTION_GENOME_ALTERED (n=365) beside it.
#
# Run from the repo root:   Rscript code/fetch_tcga_ucec_2013.R
# Idempotent; pass --force to re-download. Deps: jsonlite + tidyverse (same as the
# 2018 fetcher). The optional `curl` package turns the segment pull into ONE POST
# (~5 s); without it the script falls back to a per-sample GET loop (~2 min).
# Then knit analysis/07_tmb.Rmd (Part 3).
# =============================================================================

suppressPackageStartupMessages({
  ok <- requireNamespace("jsonlite", quietly = TRUE) && requireNamespace("readr", quietly = TRUE) &&
        requireNamespace("dplyr", quietly = TRUE) && requireNamespace("tidyr", quietly = TRUE)
  if (!ok) stop("Need jsonlite + tidyverse. Run: renv::install('jsonlite')")
  library(dplyr); library(tidyr); library(readr); library(jsonlite)
})

force   <- "--force" %in% commandArgs(TRUE)
root    <- tryCatch(here::here(), error = function(e) getwd())
api     <- "https://www.cbioportal.org/api"
study   <- "ucec_tcga_pub"                       # the 2013 Nature publication cohort
pancan  <- "ucec_tcga_pan_can_atlas_2018"        # only for ANEUPLOIDY_SCORE
dest    <- file.path(root, "data", "tcga_2013")
f_samp  <- file.path(dest, "data_clinical_sample.txt")
f_pat   <- file.path(dest, "data_clinical_patient.txt")
f_aneu  <- file.path(dest, "pancan_aneuploidy.tsv")
f_seg   <- file.path(dest, "data_cna_hg19.seg")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

get_json <- function(url) {
  d <- tryCatch(jsonlite::fromJSON(url), error = function(e)
    stop("cBioPortal API unreachable: ", conditionMessage(e)))
  if (is.null(d) || NROW(d) == 0) stop("empty response: ", url)
  as_tibble(d)
}

# clinical-data returns one LONG row per (entity, attribute); pivot to one row per entity.
fetch_clinical <- function(sid, type) {
  url <- sprintf("%s/studies/%s/clinical-data?clinicalDataType=%s&projection=SUMMARY&pageSize=10000000&pageNumber=0",
                 api, sid, type)
  message("GET ", url)
  get_json(url)
}
to_wide <- function(long, keys) {
  long |>
    distinct(across(all_of(c(keys, "clinicalAttributeId"))), .keep_all = TRUE) |>
    pivot_wider(id_cols = all_of(keys), names_from = clinicalAttributeId, values_from = value)
}

# ---------------------------------------------------------------------------
# 1-2. Clinical (sample + patient) for the 2013 study
# ---------------------------------------------------------------------------
if (all(file.exists(c(f_samp, f_pat))) && !force) {
  message("clinical already present (use --force to refresh) — skipping.")
} else {
  samp <- fetch_clinical(study, "SAMPLE")  |> to_wide(c("sampleId", "patientId")) |>
    rename(SAMPLE_ID = sampleId, PATIENT_ID = patientId)
  pat  <- fetch_clinical(study, "PATIENT") |> to_wide("patientId") |> rename(PATIENT_ID = patientId)
  write_tsv(samp, f_samp); write_tsv(pat, f_pat)
  message("Wrote ", nrow(samp), " samples x ", ncol(samp), " cols -> ", f_samp)
  message("Wrote ", nrow(pat),  " patients x ", ncol(pat), " cols -> ", f_pat)
}

# ---------------------------------------------------------------------------
# 3. ANEUPLOIDY_SCORE from the PanCancer-Atlas study (absent from the 2013 study)
# ---------------------------------------------------------------------------
if (file.exists(f_aneu) && !force) {
  message("pancan aneuploidy already present — skipping.")
} else {
  ps <- fetch_clinical(pancan, "SAMPLE") |>
    filter(clinicalAttributeId == "ANEUPLOIDY_SCORE") |>
    transmute(PATIENT_ID = patientId,
              ANEUPLOIDY_SCORE = suppressWarnings(as.numeric(value))) |>
    filter(!is.na(ANEUPLOIDY_SCORE)) |>
    distinct(PATIENT_ID, .keep_all = TRUE)
  write_tsv(ps, f_aneu)
  message("Wrote ", nrow(ps), " patient aneuploidy scores -> ", f_aneu)
}

# ---------------------------------------------------------------------------
# 4. Segmented copy number (hg19) — the vertical axis of Fig. 1a
# ---------------------------------------------------------------------------
if (file.exists(f_seg) && !force) {
  message("segments already present — skipping.")
} else {
  ids <- get_json(sprintf("%s/studies/%s/samples?pageSize=100000", api, study))$sampleId
  message("Fetching segments for ", length(ids), " samples ...")

  segs <- NULL
  # Fast path: ONE POST for the whole study (needs the `curl` package).
  if (requireNamespace("curl", quietly = TRUE)) {
    body <- jsonlite::toJSON(data.frame(studyId = study, sampleId = ids), auto_unbox = FALSE)
    h <- curl::new_handle()
    curl::handle_setopt(h, post = TRUE, postfields = body)
    curl::handle_setheaders(h, "Content-Type" = "application/json")
    segs <- tryCatch({
      r <- curl::curl_fetch_memory(paste0(api, "/copy-number-segments/fetch?projection=SUMMARY"), handle = h)
      if (r$status_code != 200) stop("HTTP ", r$status_code)
      as_tibble(jsonlite::fromJSON(rawToChar(r$content)))
    }, error = function(e) { message("  POST failed (", conditionMessage(e), ") — falling back to GET loop."); NULL })
  }
  # Fallback: one GET per sample (slower, but jsonlite-only).
  # ⚠️ pageSize IS CAPPED AT 20000 on this endpoint — a larger value returns HTTP 400
  # ("pageSize must be less than or equal to 20000"), which fromJSON() raises as an error.
  # With a tryCatch(-> NULL) that would drop the sample SILENTLY and yield a short .seg.
  # (The clinical-data endpoint has no such cap, hence the 1e7 used above — the limits
  # genuinely differ per endpoint.) Max observed is ~2000 segments/sample, so one page is
  # always enough; we still check for a full page and warn rather than truncate quietly.
  PAGE <- 20000
  if (is.null(segs)) {
    failed <- character(0); maxed <- character(0)
    segs <- lapply(seq_along(ids), function(i) {
      if (i %% 50 == 0) message("  ", i, "/", length(ids))
      r <- tryCatch(as_tibble(jsonlite::fromJSON(sprintf(
             "%s/studies/%s/samples/%s/copy-number-segments?pageSize=%d", api, study, ids[i], PAGE))),
             error = function(e) { failed <<- c(failed, ids[i]); NULL })
      if (!is.null(r) && nrow(r) >= PAGE) maxed <<- c(maxed, ids[i])
      r
    })
    segs <- bind_rows(segs[lengths(segs) > 0])
    # Never let a partial fetch masquerade as a complete one.
    if (length(failed) > 0)
      stop(length(failed), " sample(s) returned no segments (first: ", failed[[1]],
           "). Refusing to write a partial .seg — re-run, or install the `curl` package ",
           "so the single-POST path is used.")
    if (length(maxed) > 0)
      warning(length(maxed), " sample(s) hit the ", PAGE, "-row page limit — segments may be ",
              "truncated for: ", paste(utils::head(maxed, 3), collapse = ", "))
  }
  if (is.null(segs) || nrow(segs) == 0) stop("no segments returned")
  if (dplyr::n_distinct(segs$sampleId) < length(ids))
    message("note: ", length(ids) - dplyr::n_distinct(segs$sampleId), " of ", length(ids),
            " samples have no segment data (expected — not every sample was CNA-profiled).")

  # cBioPortal/Broad encode the sex chromosomes numerically (23 = X, 24 = Y); the arm
  # boundary tables use X/Y, so normalise here rather than in every downstream caller.
  out <- segs |>
    transmute(Sample       = sampleId,
              Chromosome   = dplyr::recode(as.character(chromosome), "23" = "X", "24" = "Y"),
              Start        = as.numeric(start),
              End          = as.numeric(end),
              Num_Probes   = suppressWarnings(as.integer(numberOfProbes)),
              Segment_Mean = as.numeric(segmentMean)) |>
    arrange(Sample, match(Chromosome, c(as.character(1:22), "X", "Y")), Start)
  write_tsv(out, f_seg)
  message("Wrote ", nrow(out), " segments x ", dplyr::n_distinct(out$Sample), " samples -> ", f_seg)
}

# --- verify with the project loaders -----------------------------------------
srcC <- file.path(root, "code", "attend_classes.R")
srcW <- file.path(root, "code", "load_wes_results.R")
if (file.exists(srcC) && file.exists(srcW)) {
  suppressPackageStartupMessages(source(srcC)); suppressPackageStartupMessages(source(srcW))
  cl <- load_tcga_ucec_2013()
  sg <- load_tcga_2013_seg()
  message("\nload_tcga_ucec_2013(): ", nrow(cl), " patients; ",
          "subtype called = ",   sum(!is.na(cl$subtype)), "; ",
          "CN cluster K4 = ",    sum(!is.na(cl$cn_cluster_k4)), "; ",
          "aneuploidy = ",       sum(is.finite(cl$aneuploidy)))
  message("load_tcga_2013_seg(): ", nrow(sg), " segments x ",
          dplyr::n_distinct(sg$ID), " samples")
  message("Ready — knit analysis/07_tmb.Rmd (Part 3)")
}
