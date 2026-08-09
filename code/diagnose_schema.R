#!/usr/bin/env Rscript
# =============================================================================
# diagnose_schema.R — structural correctness probe for the ATTEND pipeline
# =============================================================================
# Run on the cluster (where the real data lives) and send the WHOLE output back:
#
#     Rscript code/diagnose_schema.R 2>&1 | tee attend_schema_report.txt
#
# It prints STRUCTURE ONLY — column names, dims, category *labels* with counts,
# join overlaps, NA/range summaries, value distributions of config-critical
# columns. It deliberately does NOT print patient-level rows. (Sequencing
# barcodes appear in the join section so format mismatches can be diagnosed; if
# your governance treats barcodes as identifiers, redact that section before
# sharing — the COUNTS alone are still informative.)
#
# Each section answers a specific correctness question (see the [Q] tags). Every
# section is independently guarded, so a missing file degrades one section only.
#
# Toggles:
DO_IHC <- FALSE   # load_ihc_data() walks the whole FlowPath tree (slow); off by default.
MAF_SAMPLE_N <- 3 # how many raw .maf files to inspect for column structure.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(here); library(fs)
})

src <- function(f) tryCatch(source(here("code", paste0(f, ".R"))),
                            error = function(e) cat("!! could not source", f, ":", conditionMessage(e), "\n"))
invisible(lapply(c("load_clinical", "load_phenotypes", "load_wes_results",
                   "attend_harmonise", "attend_classes"), src))

banner <- function(txt) cat("\n", strrep("=", 78), "\n== ", txt, "\n", strrep("=", 78), "\n", sep = "")
safe   <- function(expr) tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "probe_err"))
is_err <- function(x) inherits(x, "probe_err")
show   <- function(x) print(as.data.frame(x), row.names = FALSE)   # robust for df or tibble

# ---------------------------------------------------------------------------
banner("0. Session / package availability   [Q: does the migrated stack load?]")
cat("R version:", R.version.string, "\n\n")
pkgs <- c("tidyverse","here","fs","readxl","data.table","workflowr","knitr","rmarkdown",
          "arrow","qs2","UpSetR","reclin2","naniar","skimr","pointblank",
          "survival","ggsurvfit","gtsummary","broom","ggpubr","maftools",
          "spatstat.geom","sf","dbscan","imcRtools")
inst <- rownames(installed.packages())
for (p in pkgs) cat(sprintf("  [%s] %s\n", ifelse(p %in% inst, "x", " "), p))
cat("\n([x] installed, [ ] missing — missing optional pkgs only disable their guarded sections)\n")

# ---------------------------------------------------------------------------
banner("1. Loader shapes + columns   [Q: real column names vs attend_cols placeholders]")
loaders <- list(
  clinical = load_clinical_data, imaging = load_imaging_data,
  gianlu   = load_gianlu_clinical_data,
  aneu = load_aneuploidy_scores, hrd = load_hrd_scores,
  msi  = load_msi_scores, tmb = load_tmb_scores, maf = load_maf_data
)
if (DO_IHC) loaders$ihc <- load_ihc_data

tabs <- map(loaders, ~ safe(.x()))
iwalk(tabs, ~ if (is_err(.x)) cat(sprintf("  %-9s ERROR: %s\n", .y, .x))
              else cat(sprintf("  %-9s %5d x %3d : %s\n", .y, nrow(.x), ncol(.x),
                               paste(names(.x), collapse = ", "))))

# ---------------------------------------------------------------------------
banner("2. attend_cols resolution   [Q: does every configured column exist in the data?]")
# Split each "dataset__column" config value and check the column exists in that loader.
# NB: iwalk() calls .f(value, name), so the column string is the FIRST arg.
check_col <- function(value, key) {
  if (!grepl("__", value)) return(cat(sprintf("  %-11s %-26s (no prefix — skipped)\n", key, value)))
  ds  <- sub("__.*$", "", value); col <- sub("^[^_]*__", "", value)
  tab <- tabs[[ds]]
  status <- if (is.null(tab)) "loader not run (e.g. ihc off)"
            else if (is_err(tab)) "loader ERROR"
            else if (col %in% names(tab)) "OK — found"
            else paste0("MISSING in '", ds, "' — have: ", paste(names(tab), collapse = ", "))
  cat(sprintf("  %-11s %-26s -> %s\n", key, value, status))
}
iwalk(attend_cols, check_col)

# ---------------------------------------------------------------------------
banner("3. Category vocabularies   [Q: are attend_levels / km_group tokens right?]")
g <- tabs$gianlu
cat_cols <- c("TREATMENT","MMR_STATUS","MSI_STATUS","STATUS","HISTO_GRADE","FIGO")
if (is_err(g) || is.null(g)) cat("  gianlu unavailable\n") else
  for (c in cat_cols) if (c %in% names(g)) { cat("  --", c, ":\n"); show(count(g, .data[[c]])) }
# Cross-source MSI/MMR label check (do the spellings agree across tables?)
for (c in c("MSI_STATUS","MMR_STATUS"))
  for (src_nm in c("clinical","imaging"))
    if (!is_err(tabs[[src_nm]]) && !is.null(tabs[[src_nm]])) {
      hit <- intersect(c, names(tabs[[src_nm]]))
      if (length(hit)) { cat(sprintf("  -- %s in %s:\n", c, src_nm)); show(count(tabs[[src_nm]], .data[[c]])) }
    }

# ---------------------------------------------------------------------------
banner("4. Survival encoding   [Q: is PFS_EVENT text? does recode_event() map it right?]")
if (!is_err(g) && !is.null(g)) {
  ev <- attend_cols$surv_event; tm <- attend_cols$surv_time
  ev_raw <- sub("^[^_]*__", "", ev); tm_raw <- sub("^[^_]*__", "", tm)
  if (ev_raw %in% names(g)) {
    cat("  ", ev_raw, "distinct values -> recode_event():\n")
    tibble(value = unique(as.character(g[[ev_raw]]))) |>
      mutate(recoded = recode_event(value)) |> arrange(recoded) |> show()
  } else cat("  surv_event column", ev_raw, "NOT in gianlu\n")
  if (tm_raw %in% names(g)) {
    x <- suppressWarnings(as.numeric(g[[tm_raw]]))
    cat(sprintf("\n  %s: n=%d, NA=%d, range=[%.2f, %.2f], median=%.2f (confirm UNITS = months)\n",
                tm_raw, length(x), sum(is.na(x)), min(x,na.rm=TRUE), max(x,na.rm=TRUE), median(x,na.rm=TRUE)))
  } else cat("  surv_time column", tm_raw, "NOT in gianlu\n")
}

# ---------------------------------------------------------------------------
banner("4b. aneu extra columns   [Q: which is THE score? is aneu's own survival cleaner?]")
a <- tabs$aneu
if (is_err(a) || is.null(a)) cat("  aneu unavailable\n") else {
  for (nc in intersect(c("aneuploidy","aneuploidy_score","time"), names(a))) {
    x <- suppressWarnings(as.numeric(a[[nc]]))
    cat(sprintf("  %-16s numeric: NA=%d range=[%.3f, %.3f] median=%.3f\n",
                nc, sum(is.na(x)), min(x,na.rm=TRUE), max(x,na.rm=TRUE), median(x,na.rm=TRUE)))
  }
  for (cc in intersect(c("status","arm","MMR_status"), names(a))) {
    cat("  --", cc, "(aneu):\n"); show(count(a, .data[[cc]]))
  }
}

# ---------------------------------------------------------------------------
banner("5. Barcode join coverage   [Q: do WES barcodes reach the crosswalk?]")
# The aneu loader does NO suffix stripping (unlike msi/tmb/hrd) — this is the
# section that reveals whether that breaks its mapping.
cfg <- safe(make_id_cfg()); cw <- safe(build_barcode_pid(g, cfg))
if (is_err(cw)) cat("  crosswalk build failed:", cw, "\n") else {
  cat(sprintf("  crosswalk: %d barcode->pid pairs, %d patients\n\n", nrow(cw), n_distinct(cw$pid)))
  cat("  first 5 crosswalk barcodes:", paste(head(unique(cw$barcode),5), collapse=" | "), "\n\n")
  for (nm in c("aneu","hrd","msi","tmb","maf")) {
    df <- tabs[[nm]]; if (is.null(df) || is_err(df)) { cat(sprintf("  %-5s unavailable\n", nm)); next }
    col <- cfg[[nm]]$col
    if (!col %in% names(df)) { cat(sprintf("  %-5s key col '%s' MISSING (have: %s)\n",
                                            nm, col, paste(names(df), collapse=", "))); next }
    ids <- unique(as.character(df[[col]]))
    mapped <- sum(ids %in% cw$barcode); unmapped <- setdiff(ids, cw$barcode)
    cat(sprintf("  %-5s key='%s'  %d ids, %d mapped, %d unmapped\n",
                nm, col, length(ids), mapped, length(unmapped)))
    if (length(unmapped)) cat("        e.g. unmapped:", paste(head(unmapped,5), collapse=" | "), "\n")
  }
}

# ---------------------------------------------------------------------------
banner("6. Crosswalk vs canonical attend_barcodes.csv   [Q: derived map agrees with canonical?]")
conv <- safe(read.csv(here("data","attend_barcodes.csv")))
if (is_err(conv)) cat("  attend_barcodes.csv:", conv, "\n") else {
  cat("  canonical columns:", paste(names(conv), collapse=", "), "\n")
  cand <- c("enrolment_id","enrollment_id","internal_id",
            "pid","PID","PATIENT_ID","patient_id","Patient_ID","patient","ID","Unique_Subject_Identifier")
  pidc <- setdiff(intersect(cand, names(conv)), "barcode")[1]
  cat("  detected canonical patient-id column:", ifelse(is.na(pidc), "(none)", pidc), "\n")
  if (!is_err(cw) && !is.na(pidc) && "barcode" %in% names(conv)) {
    canon <- transmute(conv, barcode = as.character(barcode), pid_canon = as.character(.data[[pidc]]))
    cmp <- inner_join(cw, canon, by = "barcode")
    cat(sprintf("  shared barcodes: %d | derived-vs-canonical pid disagreements: %d (must be 0)\n",
                nrow(cmp), sum(cmp$pid != cmp$pid_canon)))
  }
}

# ---------------------------------------------------------------------------
banner("7. Multiplicity / collapse safety   [Q: is collapse-by-mean hiding heterogeneity?]")
if (!is_err(cw)) {
  dup <- cw |> count(barcode) |> filter(n > 1)
  mul <- cw |> count(pid)     |> filter(n > 1)
  cat(sprintf("  barcodes -> >1 patient: %d (MUST be 0)\n", nrow(dup)))
  cat(sprintf("  patients -> >1 barcode: %d (legit, but these rows get averaged)\n", nrow(mul)))
}
for (nm in c("aneu","hrd","msi","tmb")) {
  df <- tabs[[nm]]
  if (is.null(df) || is_err(df) || is_err(cw)) next
  col <- cfg[[nm]]$col; if (!col %in% names(df)) next
  pid <- cw$pid[match(as.character(df[[col]]), cw$barcode)]; pid <- pid[!is.na(pid)]
  cat(sprintf("  %-5s patients with multiple rows averaged: %d\n", nm, sum(table(pid) > 1)))
}

# ---------------------------------------------------------------------------
banner("8. MAF annotation structure   [Q: which column holds pathogenicity?]")
maf_dir <- here("data","variant_annotations")
files <- safe(dir_ls(maf_dir))
if (is_err(files) || !length(files)) cat("  no MAF files at", maf_dir, "\n") else {
  cat("  ", length(files), "MAF files. Inspecting", min(MAF_SAMPLE_N, length(files)), ":\n")
  for (f in head(files, MAF_SAMPLE_N)) {
    one <- safe(data.table::fread(f, skip = "Hugo_Symbol", nrows = 5000,
                                  na.strings = c(".","","NA")))
    if (is_err(one)) { cat("   ", path_file(f), "READ ERROR:", one, "\n"); next }
    cat("\n   file:", path_file(f), "—", ncol(one), "cols\n")
    cat("   columns:", paste(names(one), collapse=", "), "\n")
    pcols <- names(one)[map_lgl(one, ~ any(str_detect(as.character(.x),
                          regex("pathogenic", ignore_case = TRUE)), na.rm = TRUE))]
    cat("   columns whose VALUES contain 'pathogenic':",
        ifelse(length(pcols), paste(pcols, collapse=", "), "(none — check the annotation column name!)"), "\n")
    # show the distinct values of the most likely clin-sig columns
    for (cc in intersect(c("CLIN_SIG","Clinical_Significance","clin_sig","IMPACT","Consequence",
                           "Variant_Classification"), names(one)))
      { cat("     ", cc, "values:", paste(head(unique(as.character(one[[cc]])),12), collapse=" | "), "\n") }
  }
  # aggregate TP53 status the pipeline actually derives
  if (!is_err(tabs$maf)) { cat("\n  load_maf_data() output TP53 status distribution:\n")
    show(tabs$maf |> count(across(matches("TP53")))) }
}

banner("DONE — paste the full output back. Sections map 1:1 to the open [Q] questions.")
