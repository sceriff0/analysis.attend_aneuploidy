# build_master.R — the join, as a function.
#
# WHY THIS IS NOT ONLY IN A REPORT
# --------------------------------
# Every report needs the master table, and until now the only thing that could produce it
# was report 10. That imposed a build order: knit 10 first, or everything else renders as
# scaffolding. Moving the join here removes that constraint entirely — get_master() reads
# the intermediate when it exists and builds it when it does not, so ANY report can be the
# first one built and `wflow_build()` with no arguments is always correct.
#
# Report 10 still exists and still owns the *narration*: it calls build_master(verbose =
# TRUE), then audits what the join cost. The arithmetic lives here so it cannot differ
# between the report that documents it and the reports that consume it.
#
# Sourcing this file is cheap; it defines functions and reads nothing.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
})

source(here::here("code", "load_clinical.R"))
source(here::here("code", "load_phenotypes.R"))
source(here::here("code", "load_wes_results.R"))
source(here::here("code", "attend_harmonise.R"))
source(here::here("code", "attend_classes.R"))
source(here::here("code", "attend_io.R"))

#' Read every raw input through its loader.
#'
#' The seam between HPC storage and the pipeline. Nothing is filtered or transformed here —
#' these are the raw shapes as delivered.
load_raw_inputs <- function() {
  list(
    clinical_data        = load_clinical_data(),
    imaging_data         = load_imaging_data(),
    ihc_data             = load_ihc_data(),
    gianlu_clinical_data = load_gianlu_clinical_data(),
    aneu_data            = load_aneuploidy_scores(),
    hrd_data             = load_hrd_scores(),
    msi_data             = load_msi_scores(),
    tmb_data             = load_tmb_scores(),
    maf_data             = load_maf_data(),
    # The canonical barcode map. Absent in a fresh checkout, so this is guarded: the join
    # itself does not need it (the crosswalk is derived from the clinical table); only
    # report 10's validation section does.
    conversion_data      = tryCatch(
      utils::read.csv(here::here("data", "attend_barcodes.csv")),
      error = function(e) NULL)
  )
}

#' The two crosswalks and the id -> pid translator, from already-loaded inputs.
build_crosswalks <- function(raw) {
  id_cfg <- make_id_cfg()
  cw_barcode_pid <- build_barcode_pid(raw$gianlu_clinical_data, id_cfg)
  cw_image_pid   <- build_image_pid(raw$imaging_data, id_cfg)
  # recover = TRUE promotes an orphaned id when its normalised key (upper-cased, trimmed)
  # maps to exactly ONE patient. Ambiguous keys are refused, never guessed — report 10's
  # recovery ledger shows precisely who is reclaimed and who is refused.
  pid_vector <- make_pid_vector(cw_barcode_pid, cw_image_pid, recover = TRUE)
  tables <- list(
    tmb = raw$tmb_data, msi = raw$msi_data, maf = raw$maf_data, hrd = raw$hrd_data,
    aneu = raw$aneu_data, gianlu = raw$gianlu_clinical_data,
    imaging = raw$imaging_data, ihc = raw$ihc_data
  )
  list(id_cfg = id_cfg, cw_barcode_pid = cw_barcode_pid, cw_image_pid = cw_image_pid,
       pid_vector = pid_vector, tables = tables,
       sets = build_sets(tables, id_cfg, pid_vector))
}

#' Build the master table: join, collapse, mutation status, TMB recompute.
#'
#' @param raw output of load_raw_inputs(); loaded here when NULL
#' @param cw output of build_crosswalks(); built here when NULL
#' @param write write the result (and each raw table) to output/clean_data/
#' @param verbose narrate, for report 10
#' @return the master tibble, one row per patient
build_master <- function(raw = NULL, cw = NULL, write = TRUE, verbose = FALSE) {
  if (is.null(raw)) raw <- load_raw_inputs()
  if (is.null(cw))  cw  <- build_crosswalks(raw)
  say <- function(...) if (verbose) cat(...)

  if (verbose) {
    spine <- names(which.max(lengths(cw$sets)))
    say("Largest single dataset:", spine, "—", length(cw$sets[[spine]]), "patients\n")
  }

  # Union of every patient any dataset reached, then left_join each collapsed dataset onto
  # it, so an absent modality becomes NA rather than dropping the patient.
  master <- tibble(pid = Reduce(union, cw$sets))
  for (nm in c("gianlu", "imaging", "ihc", "tmb", "msi", "hrd", "aneu")) {
    src <- switch(nm,
      gianlu = raw$gianlu_clinical_data, imaging = raw$imaging_data, ihc = raw$ihc_data,
      tmb = raw$tmb_data, msi = raw$msi_data, hrd = raw$hrd_data, aneu = raw$aneu_data)
    master <- left_join(master,
      collapse_pid(src, cw$id_cfg[[nm]], nm, cw$pid_vector), by = "pid")
  }

  # The MAF is deliberately excluded from the averaging: the mean of a genomic position is
  # not a genomic position. Mutation status is derived per patient as booleans instead.
  master <- master |>
    left_join(pathogenic_by_patient(raw$maf_data, cw$cw_barcode_pid), by = "pid")

  for (s in names(cw$sets)) master[[paste0("in_", s)]] <- master$pid %in% cw$sets[[s]]

  # Ancestry-aware TMB recompute, anchored to the upstream TMB_SCORE scale. Knit-safe: if
  # the gnomAD-annotated MAFs are not synced, load_maf_tmb() returns empty and the
  # recomputed columns are simply absent.
  maf_tmb <- tryCatch(load_maf_tmb(), error = function(e) tibble())
  if (nrow(maf_tmb) > 0) {
    anchor <- if (attend_cols$tmb %in% names(master))
      master |> transmute(pid, tmb_anchor = suppressWarnings(as.numeric(.data[[attend_cols$tmb]]))) else NULL
    tmb_vars <- tmb_variant_table(maf_tmb, cw$cw_barcode_pid,
                                  anc_by_pid = ancestry_by_pid(master),
                                  anchor_tmb = anchor)
    if (!is.null(tmb_vars)) master <- left_join(master, tmb_vars, by = "pid")
  }
  master <- apply_nassar_recalibration(master, ancestry_by_pid(master))

  if (verbose) {
    present <- tmb_defs_present(master)
    say("TMB definitions available downstream:\n")
    print(tibble(definition = tmb_def_labels(present), column = unname(present)))
  }

  if (write) {
    p <- write_intermediate(master, "attend_master_joined")
    say("Wrote", basename(p), ":", nrow(master), "patients x", ncol(master), "columns\n")
    # Each raw table alongside the master, so the audit can reconstruct the join without
    # re-reading HPC storage.
    datasets <- raw[!vapply(raw, is.null, logical(1))]
    if (nrow(maf_tmb) > 0) datasets$maf_tmb <- maf_tmb
    iwalk(datasets, ~ write_intermediate(.x, .y))
  }
  master
}

#' The master table, built on demand.
#'
#' THIS is what removes the build order. Reads the intermediate when it exists; builds and
#' writes it when it does not. Every report calls this instead of read_intermediate(), so
#' no report depends on another having been knitted first.
#'
#' @param refresh force a rebuild even when the intermediate exists
get_master <- function(refresh = FALSE) {
  if (!refresh) {
    m <- tryCatch(read_intermediate("attend_master_joined"), error = function(e) NULL)
    if (!is.null(m) && nrow(m) > 0) return(m)
  }
  message("master not found — building it now (code/build_master.R)")
  build_master(write = TRUE, verbose = FALSE)
}

# Allow `Rscript code/build_master.R` to rebuild it from the shell.
if (sys.nframe() == 0L) invisible(build_master(write = TRUE, verbose = TRUE))
