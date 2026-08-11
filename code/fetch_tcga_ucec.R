#!/usr/bin/env Rscript
# =============================================================================
# fetch_tcga_ucec.R — R-native fetch of the TCGA UCEC clinical table (legacy
# PanCanAtlas export; the 2013 publication cohort fetched by fetch_tcga_ucec_2013.R
# is what report 41 actually uses)
#
# Cohort B (PRIMARY tumours) for the ATTEND-vs-TCGA MMRd aneuploidy contrast. Pulls ONLY
# the clinical attributes needed from the cBioPortal REST API — a few hundred KB,
# NOT the ~200 MB study tarball — pivots them to wide, and writes the two files
# load_tcga_ucec() reads into data/tcga/. Idempotent; pass --force to re-download.
#
# Run from the repo root:   Rscript code/fetch_tcga_ucec.R
# Deps: tidyverse + jsonlite (jsonlite is the only extra; renv::install("jsonlite")).
# load_tcga_ucec() is the fallback TCGA source; report 41 prefers load_tcga_ucec_2013().
#
# Attributes fetched (confirmed present for this study, 2026-07):
#   SAMPLE-level : ANEUPLOIDY_SCORE, MSI_SCORE_MANTIS, MSI_SENSOR_SCORE,
#                  FRACTION_GENOME_ALTERED, TMB_NONSYNONYMOUS
#   PATIENT-level: SUBTYPE, OS_MONTHS, OS_STATUS   (SUBTYPE = UCEC_MSI/_POLE/_CN_*)
# (The API returns every attribute; the loader just reads the ones it needs by name.)
# =============================================================================

suppressPackageStartupMessages({
  ok <- requireNamespace("jsonlite", quietly = TRUE) && requireNamespace("readr", quietly = TRUE) &&
        requireNamespace("dplyr", quietly = TRUE) && requireNamespace("tidyr", quietly = TRUE)
  if (!ok) stop("Need jsonlite + tidyverse. Run: renv::install('jsonlite')  (tidyverse should already be present)")
  library(dplyr); library(tidyr); library(readr); library(jsonlite)
})

force  <- "--force" %in% commandArgs(TRUE)
root   <- tryCatch(here::here(), error = function(e) getwd())
study  <- "ucec_tcga_pan_can_atlas_2018"
api    <- "https://www.cbioportal.org/api"
dest   <- file.path(root, "data", "tcga")
f_samp <- file.path(dest, "data_clinical_sample.txt")
f_pat  <- file.path(dest, "data_clinical_patient.txt")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

if (all(file.exists(c(f_samp, f_pat))) && !force) {
  message("data/tcga/ already populated (use --force to refresh). Skipping download.")
} else {
  fetch_clinical <- function(type) {                       # type = "SAMPLE" | "PATIENT"
    url <- sprintf("%s/studies/%s/clinical-data?clinicalDataType=%s&projection=SUMMARY&pageSize=10000000&pageNumber=0",
                   api, study, type)
    message("GET ", url)
    dat <- tryCatch(jsonlite::fromJSON(url), error = function(e) stop("cBioPortal API unreachable: ", conditionMessage(e)))
    if (is.null(dat) || NROW(dat) == 0) stop("empty response for ", type)
    as_tibble(dat)
  }
  to_wide <- function(long, keys) {                        # long -> one row per key, one col per attribute
    long |>
      distinct(across(all_of(c(keys, "clinicalAttributeId"))), .keep_all = TRUE) |>
      pivot_wider(id_cols = all_of(keys), names_from = clinicalAttributeId, values_from = value)
  }

  samp <- fetch_clinical("SAMPLE")  |> to_wide(c("sampleId", "patientId")) |>
    rename(SAMPLE_ID = sampleId, PATIENT_ID = patientId)
  pat  <- fetch_clinical("PATIENT") |> to_wide("patientId") |> rename(PATIENT_ID = patientId)

  write_tsv(samp, f_samp)
  write_tsv(pat,  f_pat)
  message("Wrote ", nrow(samp), " samples x ", ncol(samp), " cols -> ", f_samp)
  message("Wrote ", nrow(pat),  " patients x ", ncol(pat), " cols -> ", f_pat)
}

# --- verify with the project loader ------------------------------------------
srcC <- file.path(root, "code", "attend_classes.R")
srcW <- file.path(root, "code", "load_wes_results.R")
if (file.exists(srcC) && file.exists(srcW)) {
  suppressPackageStartupMessages(source(srcC)); suppressPackageStartupMessages(source(srcW))
  t <- load_tcga_ucec()
  message("\nload_tcga_ucec(): ", nrow(t), " samples; MMRd = ", sum(t$mmrd, na.rm = TRUE),
          "; aneuploidy non-NA = ", sum(is.finite(t$aneuploidy)),
          "; with OS = ", sum(!is.na(t$os_months)))
  message("Ready — load_tcga_ucec() can now read data/tcga_ucec/")
}

# --- Fallbacks (if the REST API is down) -------------------------------------
# (a) Bulk tarball:  ./code/fetch_tcga_ucec.sh   (streams the study, extracts clinical files)
# (b) Bioconductor:  BiocManager::install("cBioPortalData")
#     mae <- cBioPortalData::cBioDataPack("ucec_tcga_pan_can_atlas_2018")
#     readr::write_tsv(as.data.frame(MultiAssayExperiment::colData(mae)), "data/tcga/data_clinical_sample.txt")
