# =============================================================================
# attend_classes.R
# Shared configuration + class-derivation helpers for the ATTEND analyses.
#
# Sourced by:
#   analysis/01-data-integration.Rmd       (mutation status for the master)
#   analysis/04-survival.Rmd
#   analysis/05-aneuploidy-mmrd.Rmd
#
# >>> SET YOUR COLUMN NAMES IN `attend_cols` BELOW. <<<
# Everything downstream reads from here, so you only edit names in one place.
# Columns refer to the wide master table (output/clean_data/attend_master_joined),
# whose columns are PREFIXED by dataset (gianlu__, aneu__, tmb__, ...).
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))

# --- 1. Column names (placeholders — edit these) ----------------------------
# Verified against the real data via code/diagnose_schema.R. Names are the
# master's prefixed columns (dataset__column).
attend_cols <- list(
  pid        = "pid",
  surv_time  = "gianlu__PFS_MONTHS",       # months, range 0.03-53, no NA (verified)
  surv_event = "gianlu__PFS_EVENT",        # TEXT "No event"/"PD"/"Death" -> recode_event() => 0/1/1 (verified)
  km_group   = "gianlu__TREATMENT",        # treatment arm Atezolizumab/Placebo = the red/cyan split (verified)
  aneuploidy = "aneu__aneuploidy_score",   # continuous [0,1]; the bare `aneuploidy` column is empty (verified)
  tmb        = "tmb__TMB_SCORE",           # mut/Mb (verified)
  hrd        = "hrd__HRD_Score",           # capital S; hrd table also has LOH/TAI/LST sub-scores (verified)
  msi_status = "gianlu__MSI_STATUS",       # categorical; WES gives the continuous msi__MSI_SCORE (verified)
  mmr_status = "gianlu__MMR_STATUS"        # categorical "Deficient"/"Intact" (verified)
)
# Alternative survival source: the aneu table carries its OWN clean numeric survival
# (aneu__time + aneu__status as 0/1 + aneu__arm), keyed by barcode with 100% mapping.
# Repoint surv_time/surv_event/km_group there if you prefer numeric status over PFS text.

# --- 1b. MAF / mutation-status configuration --------------------------------
# The MAF is NOT mean-collapsed into the master (averaging variant rows or
# per-gene flags is meaningless). Mutation status is derived per patient by
# mutation_status_long() / pathogenic_by_patient() below, which support TWO
# shapes — set `format` or leave "auto" to detect:
#
#   "wide"  one row per sample, one column per gene's status. CURRENT ATTEND
#           data: `ID` + `TP53_status` (extensible: add `KRAS_status`, ...).
#           A gene is altered for a patient if ANY of their samples is positive
#           (see attend_mut_positive). Column -> gene by stripping `status_suffix`.
#
#   "long"  one row per variant: `ID` (barcode) + `Hugo_Symbol` +
#           `Variant_Classification`. The real MAF shape — load it through
#           maftools (read_maf_maftools() below) for proper variant analysis.
attend_maf <- list(
  format        = "auto",                  # "auto" | "wide" | "long"
  sample_col    = "ID",                    # barcode column (both shapes)
  gene_col      = "Hugo_Symbol",           # long: gene symbol column
  class_col     = "Variant_Classification",# long: consequence column
  status_suffix = "_status",               # wide: per-gene status cols are <GENE>_status
  # Protein-change column, by candidate name — annotators disagree (VEP emits
  # HGVSp_Short, ANNOVAR aaChange, older maftools Protein_Change), so the loaders
  # and every consumer resolve it with intersect() rather than assuming one.
  # Lives here, not under attend_pole, because TWO things need it now: the POLE
  # hotspot call and maftools::oncodrive(), which cannot run without it at all.
  protein_col   = c("HGVSp_Short", "Protein_Change", "amino_acid_change", "aaChange")
)

# Auto-detect of the "altered" value in wide per-gene status columns. Anything
# matching (case-insensitive, trimmed) is TRUE; NA/empty stays NA; all else FALSE.
# Includes "abnormal" because load_maf_data() emits "ABNORMAL"/"NORMAL" for TP53.
attend_mut_positive <- c(
  "mutated", "mutant", "mut", "pathogenic", "likely_pathogenic", "altered",
  "alteration", "abnormal", "positive", "pos", "yes", "y", "true", "t", "1",
  "deleterious", "loss", "loh"
)

# --- 2. Thresholds & category labels ----------------------------------------
attend_thresholds <- list(
  tmb        = 10,    # TMB-high if >= 10 mut/Mb
  hrd        = 50,    # HRD-high if >= 50
  aneuploidy = 0.1,   # aneuploidy-high if aneu__aneuploidy_score >= 0.1 (fixed cutoff; was median split)
  response_months = 6 # responder if progression-free beyond 6 months (see add_response_class)
)

attend_levels <- list(
  msi_unstable  = "Instable (MSI)",  # gianlu MSI_STATUS unstable label (verified; stable = "Stable (MSS)")
  mmr_deficient = "Deficient"        # gianlu MMR_STATUS deficient label (verified; proficient = "Intact")
)

# --- Manual patient highlight (every per-patient point/box plot) -------------------
# A place to flag SPECIFIC patients of interest with a SPECIFIC colour on top of any
# per-patient figure. Each entry is a named subset: `ids` — patient `pid`s, sequencing
# barcodes, OR image ids, all three resolved through the crosswalks (see below) — drawn on
# TOP of the base layer in `color`, larger, so the subset stands out from the fill.
# aneuploidy-HIGH tumours are ALREADY marked by point SHAPE (triangle); this is an
# orthogonal manual override. Add groups like:
#   attend_highlight$groups$responders <- list(ids = c("P07","P12"), color = "#009988")
#
# ONE GROUP: `polipo`, the single sample 21S188, drawn in the palette's reserved blue.
# It is highlighted on EVERY report that draws per-patient points, with NO per-report
# opt-out. The point of marking one sample is that a reader can locate it in every
# distribution it appears in; a figure that silently drops it is worse than one that never
# carried it, because the reader cannot tell the difference. Reports 05 and 06 previously
# wrapped highlight_points() to switch the overlay off — that is gone, deliberately.
#
# THERE IS NO `cohort` GROUP, and it must not come back as one. It marked the patients
# holding IHC / imaging data in yellow, and it failed on both counts a highlight has to
# meet. It covered roughly 16 of ~40 patients, so "highlighted" stopped meaning "look at
# this one" and became a second competing fill fighting the semantic palette underneath it;
# and #F0E442 yellow was the least visible overlay in the set on exactly the pale boxes it
# landed on most often. It was also not a biological group at all — it is data
# availability, which the master already carries as the `in_*` membership columns and which
# report 01 already draws as a set diagram. That is where a reader should learn it.
#
# 21S188 IS A SEQUENCING BARCODE, NOT A pid — AND THAT USED TO BREAK SILENTLY.
# Most figures are built from the master, which is keyed on `pid` and carries no barcode
# column, so the literal-string match found nothing and the overlay drew nothing on the
# majority of plots without erroring. The ids are therefore expanded through the SAME
# crosswalks the join uses: build_crosswalks() and get_master() both call
# register_highlight_xwalk() (code/attend_plots.R), and highlight_group_of() matches on the
# expanded set. Configure an id in whichever space you know it in; the resolution is the
# pipeline's problem, not the analyst's.
attend_highlight <- list(
  groups = list(
    # attend_plots.R rule [C] reserves `blue` (#0077BB) for exactly this. Highlight points
    # are drawn ON TOP of boxes filled with the MMR, aneuploidy, TP53 and response colours,
    # so a highlight sharing a hex is invisible on the one box it exists to mark. That is
    # what happened when polipo was #CC79A7 (identical to non-responder, so the point
    # vanished on report 06's non-responder box) and when the retired cohort group was
    # #E69F00 (identical to MMRd, same problem on report 03). test_figure_system.R and
    # test_plot_style.R rule [5] both pin it.
    polipo = list(ids = c("21S188"), color = "#0077BB")
  )
)

# ONE switch, pipeline-wide, for whether the highlight MARK is drawn.
#
# FALSE draws every patient in the plain base colour and omits the blue polipo point and
# its legend. The base points themselves are unaffected — at the call sites that pass
# `points = FALSE` to attend_box(), highlight_points() IS the points layer, so switching
# the mark off must not take the cohort with it.
#
# It lives here, as a single global, and NOT as an argument a report passes. The ban
# test_figure_system.R enforces is on a report opting out on its own: reports 05 and 06
# each used to carry a wrapper that dropped the overlay, so the mark was present in some
# distributions and absent in others, and a reader could not tell "not this patient" from
# "not this figure". A global keeps every report saying the same thing, whichever way it
# is set. Flip it here and the whole site follows.
attend_highlight_show <- FALSE

# --- 3. Gene panels ---------------------------------------------------------
# `attend_gene_panel` is THE MMR PANEL and nothing else. It is used two ways:
# one gene at a time (per-gene <GENE>-mut/<GENE>-wt) and as a union (>=1
# pathogenic variant in ANY of these => "panel-altered"). The union also feeds
# `panel_pathogenic` on the master, so the panel's membership defines what that
# column means — keep it to genes whose loss actually causes mismatch repair
# deficiency. MSH3 is a canonical MMR gene (MutSbeta, with MSH2) and belongs here.
attend_gene_panel <- c("PMS2", "MLH1", "MSH2", "MSH6", "MSH3")

# Wnt / beta-catenin pathway genes. Reported ALONGSIDE the MMR panel gene-by-gene,
# but deliberately kept OUT of `attend_gene_panel`: CTNNB1 and APC are not mismatch
# repair genes, so folding them into the union would silently redefine
# `panel_pathogenic` (and every "panel-altered" comparison built on it) to mean
# "MMR or Wnt altered". They get their own per-gene panels; there is no Wnt union.
attend_gene_panel_wnt <- c("CTNNB1", "APC")

# Every gene plotted ONE AT A TIME in the mutation-vs-aneuploidy report, in display
# order: MMR panel first, then Wnt. TP53 is handled separately (it is already a
# first-class column on the master, `TP53_pathogenic`).
attend_genes_per_gene <- c(attend_gene_panel, attend_gene_panel_wnt)

# Variant_Classification values treated as pathogenic for the LONG MAF only.
# NOTE: Missense is included here but usually needs a real pathogenicity annotation
# (ClinVar / OncoKB). Drop it or switch to an annotation column if you have one.
attend_pathogenic_classes <- c(
  "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation",
  "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation",
  "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"
)

# --- 4. Small guard ---------------------------------------------------------
# Returns TRUE only if all named columns exist; otherwise messages which are missing.
have_cols <- function(df, ...) {
  need <- c(...)
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    message("Missing column(s) — set them in attend_cols: ", paste(miss, collapse = ", "))
    return(FALSE)
  }
  TRUE
}

# --- 5. Derive molecular classes on a patient-level data frame --------------
# Adds: TMB_class, HRD_class, aneuploidy_class, MSI_class, MMR_class.
# Missing source columns yield all-NA classes (so it never errors before you set names).
add_molecular_classes <- function(df,
                                  cols = attend_cols,
                                  thr  = attend_thresholds,
                                  lev  = attend_levels) {
  n   <- nrow(df)
  num <- function(col) if (col %in% names(df)) suppressWarnings(as.numeric(df[[col]])) else rep(NA_real_, n)
  chr <- function(col) if (col %in% names(df)) as.character(df[[col]]) else rep(NA_character_, n)

  a       <- num(cols$aneuploidy)
  a_cut   <- if (is.null(thr$aneuploidy)) stats::median(a, na.rm = TRUE) else thr$aneuploidy

  df$TMB_class        <- ifelse(num(cols$tmb) >= thr$tmb, "TMB-high", "TMB-low")
  df$HRD_class        <- ifelse(num(cols$hrd) >= thr$hrd, "HRD-high", "HRD-low")
  df$aneuploidy_class <- ifelse(a >= a_cut, "aneuploidy-high", "aneuploidy-low")
  df$MSI_class        <- ifelse(chr(cols$msi_status) == lev$msi_unstable, "MSI-high", "MSS")
  df$MMR_class        <- chr(cols$mmr_status)
  df
}

# Recode a free-text survival event to the 1=event / 0=censored Surv() needs.
# The clinical PFS_EVENT is text ("PD"/"Death"/"no event"); these tokens -> 1, the
# rest -> 0, unknown/empty -> NA. Idempotent for data already coded 0/1 or T/F, so
# the survival helpers apply it unconditionally.
attend_event_positive <- c(
  "pd", "death", "dead", "deceased", "progression", "progressed",
  "relapse", "recurrence", "event", "1", "true", "yes"
)
recode_event <- function(x, positive = attend_event_positive) {
  v   <- tolower(trimws(as.character(x)))
  out <- as.integer(v %in% tolower(positive))
  out[is.na(x) | v %in% c("", "na", "nan", "none")] <- NA_integer_
  out
}

# --- 5b. Treatment response, as a landmark on the survival endpoint ---------
# "Responder" is defined on time alone: still progression-free past
# thr$response_months (attend_cols$surv_time / $surv_event, i.e. PFS by default).
#
# The third category is the one that matters. A patient who progressed at or before the
# landmark is a NON-RESPONDER; a patient who passed it event-free is a RESPONDER; but a
# patient whose follow-up merely STOPPED before the landmark (censored, no event) has no
# observable response status at all. Coding those as non-responders counts lost follow-up
# as treatment failure and biases the responder arm upward -- the standard failure of a
# naive landmark split. They are returned as NA so a report can print how many were lost
# rather than absorb them into a group.
#
# An event AFTER the landmark does not demote a responder: the question is what happened
# by month `response_months`, not eventual outcome. Adds response_time / response_event /
# response_class; missing source columns give an all-NA class, like add_molecular_classes().
add_response_class <- function(df, cols = attend_cols, thr = attend_thresholds) {
  n     <- nrow(df)
  time  <- if (cols$surv_time  %in% names(df)) suppressWarnings(as.numeric(df[[cols$surv_time]]))
           else rep(NA_real_, n)
  event <- if (cols$surv_event %in% names(df)) recode_event(df[[cols$surv_event]])
           else rep(NA_integer_, n)
  cut   <- thr$response_months

  cls <- rep(NA_character_, n)
  cls[!is.na(time) & time >  cut] <- "responder"
  cls[!is.na(time) & time <= cut & !is.na(event) & event == 1L] <- "non-responder"

  df$response_time  <- time
  df$response_event <- event
  df$response_class <- factor(cls, levels = c("responder", "non-responder"))
  df
}

# Why each patient did or did not get a response_class, as a one-row-per-reason table.
# Reported next to every response figure so the denominator is never implicit: the
# censored-early count is the cohort this definition cannot speak about.
response_counts <- function(df, thr = attend_thresholds) {
  cut <- thr$response_months
  tibble::tibble(
    status = c(paste0("responder (progression-free > ", cut, " mo)"),
               paste0("non-responder (event <= ", cut, " mo)"),
               paste0("unclassifiable (censored <= ", cut, " mo)"),
               "unclassifiable (no survival time)"),
    n = c(sum(df$response_class == "responder",     na.rm = TRUE),
          sum(df$response_class == "non-responder", na.rm = TRUE),
          sum(is.na(df$response_class) & !is.na(df$response_time)),
          sum(is.na(df$response_time)))
  )
}

# --- 6. Mutation-derived classes from the MAF -------------------------------
# Format-aware: works for the wide per-gene-status file (current ATTEND data) and
# for a long per-variant MAF. Needs a barcode->pid crosswalk (build it with
# build_barcode_pid() from attend_harmonise.R). Keeping this OUT of the master's
# mean-collapse path is deliberate — averaging mutation flags is meaningless.

# TRUE / FALSE / NA test for a wide status value, against the positive-token set.
is_mut_positive <- function(x, positive = attend_mut_positive) {
  v   <- tolower(trimws(as.character(x)))
  out <- v %in% tolower(positive)
  out[is.na(x) | v %in% c("", "na", "nan", "none")] <- NA
  out
}

# Decide whether a MAF data frame is "wide" or "long".
detect_maf_format <- function(maf, maf_cfg = attend_maf) {
  if (!identical(maf_cfg$format, "auto")) return(maf_cfg$format)
  if (all(c(maf_cfg$gene_col, maf_cfg$class_col) %in% names(maf))) return("long")
  "wide"   # default: a sample column + per-gene status columns
}

# Tidy per-patient, per-gene mutation status: one row per (pid, gene) with a
# logical `mutated` (TRUE if ANY of the patient's samples is altered for that gene).
mutation_status_long <- function(maf, crosswalk,
                                 cols               = attend_cols,
                                 maf_cfg            = attend_maf,
                                 pathogenic_classes = attend_pathogenic_classes,
                                 positive           = attend_mut_positive) {
  fmt <- detect_maf_format(maf, maf_cfg)

  if (fmt == "long") {
    maf |>
      transmute(barcode = as.character(.data[[maf_cfg$sample_col]]),
                gene    = as.character(.data[[maf_cfg$gene_col]]),
                vclass  = as.character(.data[[maf_cfg$class_col]])) |>
      inner_join(crosswalk, by = "barcode") |>
      filter(!is.na(pid)) |>
      group_by(pid, gene) |>
      summarise(mutated = any(vclass %in% pathogenic_classes), .groups = "drop")
  } else {
    sample_col  <- maf_cfg$sample_col
    status_cols <- grep(paste0(maf_cfg$status_suffix, "$"), names(maf), value = TRUE)
    if (length(status_cols) == 0)               # no `_status` suffix: take every payload col
      status_cols <- setdiff(names(maf), sample_col)

    maf |>
      transmute(barcode = as.character(.data[[sample_col]]), across(all_of(status_cols))) |>
      pivot_longer(all_of(status_cols), names_to = "gene", values_to = "status") |>
      mutate(gene    = sub(paste0(maf_cfg$status_suffix, "$"), "", gene),
             mutated = is_mut_positive(status, positive)) |>
      inner_join(crosswalk, by = "barcode") |>
      filter(!is.na(pid)) |>
      group_by(pid, gene) |>
      summarise(mutated = any(mutated, na.rm = TRUE), .groups = "drop")
  }
}

# One row per patient: pid, TP53_pathogenic, panel_pathogenic (logical).
# Patients in the crosswalk but with no altered gene are FALSE (not dropped).
pathogenic_by_patient <- function(maf, crosswalk,
                                  gene_panel         = attend_gene_panel,
                                  cols               = attend_cols,
                                  maf_cfg            = attend_maf,
                                  pathogenic_classes = attend_pathogenic_classes,
                                  positive           = attend_mut_positive) {
  ml <- mutation_status_long(maf, crosswalk, cols, maf_cfg, pathogenic_classes, positive)

  tp53  <- ml |> filter(gene == "TP53",       mutated) |> distinct(pid) |> mutate(TP53_pathogenic  = TRUE)
  panel <- ml |> filter(gene %in% gene_panel, mutated) |> distinct(pid) |> mutate(panel_pathogenic = TRUE)

  crosswalk |> distinct(pid) |>
    left_join(tp53,  by = "pid") |>
    left_join(panel, by = "pid") |>
    mutate(TP53_pathogenic  = replace_na(TP53_pathogenic,  FALSE),
           panel_pathogenic = replace_na(panel_pathogenic, FALSE))
}

# Optional: read a REAL per-variant MAF with maftools and return it in the long
# shape mutation_status_long() expects (ID / Hugo_Symbol / Variant_Classification).
# maftools::read.maf drops silent variants by default — fine, we only count
# pathogenic consequences. Install with: renv::install("bioc::maftools").
read_maf_maftools <- function(maf_path, maf_cfg = attend_maf, ...) {
  if (!requireNamespace("maftools", quietly = TRUE))
    stop("maftools not installed — renv::install('bioc::maftools'); renv::snapshot()")
  m <- maftools::read.maf(maf = maf_path, ...)
  tibble::as_tibble(m@data) |>
    transmute(!!maf_cfg$sample_col := as.character(Tumor_Sample_Barcode),
              !!maf_cfg$gene_col   := as.character(Hugo_Symbol),
              !!maf_cfg$class_col  := as.character(Variant_Classification))
}

# --- 7. Survival helpers ----------------------------------------------------
# Engine: survival. Faceted KM: survminer::ggsurvplot_facet — matches the reference
# attend_wes_analysis figures (arm strata, faceted by molecular axes, per-facet
# log-rank p). Single-curve KM: ggsurvfit (km_plot, optional). Pub table: gtsummary
# (surv_table). Cox tidy: broom (cox_table). The calling Rmd loads survminer.

# Single comparison KM: coloured curves for `group`, log-rank p, risk table.
# Same palette fix as km_facet() below: c("red", "turquoise3") named no palette in this project.
km_plot <- function(df, group,
                    time    = attend_cols$surv_time,
                    event   = attend_cols$surv_event,
                    palette = if (exists("attend_pal")) attend_pal(2) else c("#0077BB", "#CC3311"),
                    title   = NULL) {
  d <- df |> mutate(across(all_of(event), recode_event)) |>
    filter(if_all(all_of(c(time, event, group)), ~ !is.na(.)))
  f <- stats::as.formula(sprintf("survival::Surv(`%s`, `%s`) ~ `%s`", time, event, group))
  ggsurvfit::survfit2(f, data = d) |>
    ggsurvfit::ggsurvfit(linewidth = 1) +
    ggsurvfit::add_pvalue() +
    ggsurvfit::add_risktable() +
    ggplot2::scale_colour_manual(values = palette) +
    # ggsurvfit's default x label is a bare "time". Survival is in MONTHS throughout this
    # pipeline (attend_thresholds$response_months; surv_table(times = c(12, 24, 36, 60))).
    ggplot2::labs(title = title, x = "Time (months)", y = "Survival probability")
}

# Faceted KM matching the reference attend_wes_analysis figures: survminer's
# ggsurvplot_facet draws `group` as coloured curves (treatment arm by default),
# `facets` (1-2 vars, e.g. c("aneuploidy_class","MMR_class")) split the grid, and
# each panel carries its own log-rank p (+ test name) at fixed coordinates.
# The event column is recoded to 0/1 first (PFS_EVENT is free text).
# `palette` was c("red", "turquoise3") — two base-R colour names that belong to no palette
# in this project, so every KM curve in report 04 was the one figure family that could not
# be read against the rest of the site. It now takes the project palette (attend_plots.R),
# with a fallback so attend_classes.R stays sourceable on its own.
km_facet <- function(df, group, facets,
                     time    = attend_cols$surv_time,
                     event   = attend_cols$surv_event,
                     palette = if (exists("attend_pal")) attend_pal(2) else c("#0077BB", "#CC3311"),
                     title   = NULL) {
  keep_cols <- c(time, event, group, facets)
  d <- df |> mutate(across(all_of(event), recode_event)) |>
    filter(if_all(all_of(keep_cols), ~ !is.na(.)))
  if (nrow(d) == 0) { message("km_facet(): no complete rows for ", title); return(invisible(NULL)) }
  # survminer indexes data[, facet.by]; on a tibble a SINGLE facet returns a 1-col
  # tibble (-> "cannot xtfrm data frames"). A base data.frame returns a vector.
  d <- as.data.frame(d)

  f    <- stats::as.formula(sprintf("survival::Surv(`%s`, `%s`) ~ `%s`", time, event, group))
  fit  <- survival::survfit(f, data = d)
  fit$call$formula <- f   # survminer reads fit$call$formula — must be the evaluated formula, not `f`
  tmax <- max(suppressWarnings(as.numeric(d[[time]])), na.rm = TRUE)

  survminer::ggsurvplot_facet(
    fit, data = d, facet.by = facets,
    palette           = palette,
    pval              = TRUE, pval.method = TRUE,
    pval.coord        = c(tmax * 0.5, 0.95),
    pval.method.coord = c(tmax * 0.5, 0.80),
    censor.shape      = "+",
    short.panel.labs  = TRUE,
    # Survival time is in MONTHS everywhere in this pipeline (attend_thresholds$response_months,
    # surv_table(times = c(12, 24, 36, 60))). The axis said only "Time".
    xlab = "Time (months)", ylab = "Survival probability", legend.title = group
  ) + ggplot2::labs(title = title)
}

# Publication survival table (gtsummary): survival at the given times by `group`.
surv_table <- function(df, group, times = c(12, 24, 36, 60),
                       time  = attend_cols$surv_time,
                       event = attend_cols$surv_event) {
  d <- df |> mutate(across(all_of(event), recode_event)) |>
    filter(if_all(all_of(c(time, event, group)), ~ !is.na(.)))
  f   <- stats::as.formula(sprintf("survival::Surv(`%s`, `%s`) ~ `%s`", time, event, group))
  fit <- survival::survfit(f, data = d)
  fit$call$formula <- f   # gtsummary/cardx requires an *evaluated* formula in the survfit call
  gtsummary::tbl_survfit(fit, times = times, label_header = "**{time}-month survival**")
}

# Cox model summary: tidy hazard ratios (broom) + a gtsummary regression table.
cox_table <- function(df, covariates,
                      time  = attend_cols$surv_time,
                      event = attend_cols$surv_event) {
  rhs <- paste(sprintf("`%s`", covariates), collapse = " + ")
  f   <- stats::as.formula(sprintf("survival::Surv(`%s`, `%s`) ~ %s", time, event, rhs))
  d   <- df |> mutate(across(all_of(event), recode_event))
  fit <- survival::coxph(f, data = d)
  list(tidy  = broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE),
       table = gtsummary::tbl_regression(fit, exponentiate = TRUE))
}

# Log-rank + single-covariate Cox effect of `group` on survival, fully guarded.
# Returns a ONE-ROW tibble: n, n_events, n_groups, logrank_p, HR, HR_lo, HR_hi
# (HR = second level vs the reference first level of `group`). Every field is NA
# when the comparison is undefined — `group` has <2 levels present, zero events,
# or fewer than `min_n` complete rows — so a threshold sweep degrades to NA at the
# ends instead of erroring. Event text is recoded 0/1 first (recode_event), exactly
# as km_facet()/km_plot() do. This is the shared primitive behind any cutpoint-sweep
# sensitivity sweeps: the per-threshold KM GRID and the SUMMARY curve read their
# numbers from the same code path and therefore cannot drift.
logrank_cox_effect <- function(df, group,
                               time  = attend_cols$surv_time,
                               event = attend_cols$surv_event,
                               min_n = 5) {
  na_row <- tibble::tibble(n = 0L, n_events = 0L, n_groups = 0L,
                           logrank_p = NA_real_, HR = NA_real_,
                           HR_lo = NA_real_, HR_hi = NA_real_)
  if (!all(c(time, event, group) %in% names(df))) return(na_row)
  if (!requireNamespace("survival", quietly = TRUE)) return(na_row)

  d <- df |>
    dplyr::mutate(dplyr::across(dplyr::all_of(event), recode_event)) |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(time, event, group)), ~ !is.na(.)))
  d[[group]] <- droplevels(as.factor(d[[group]]))

  n <- nrow(d); n_ev <- sum(d[[event]] == 1, na.rm = TRUE)
  n_grp <- nlevels(d[[group]])
  if (n < min_n || n_ev < 1L || n_grp < 2L)
    return(tibble::tibble(n = n, n_events = n_ev, n_groups = n_grp,
                          logrank_p = NA_real_, HR = NA_real_,
                          HR_lo = NA_real_, HR_hi = NA_real_))

  f  <- stats::as.formula(sprintf("survival::Surv(`%s`, `%s`) ~ `%s`", time, event, group))
  lr <- tryCatch(survival::survdiff(f, data = d), error = function(e) NULL)
  logrank_p <- if (is.null(lr)) NA_real_
               else stats::pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)

  # Cox HR only when the group is binary (a single interpretable HR); >2 levels -> NA HR.
  hr <- hr_lo <- hr_hi <- NA_real_
  if (n_grp == 2L) {
    cx <- tryCatch(summary(survival::coxph(f, data = d)), error = function(e) NULL)
    if (!is.null(cx)) {
      hr    <- unname(cx$conf.int[1, "exp(coef)"])
      hr_lo <- unname(cx$conf.int[1, "lower .95"])
      hr_hi <- unname(cx$conf.int[1, "upper .95"])
    }
  }
  tibble::tibble(n = n, n_events = n_ev, n_groups = n_grp,
                 logrank_p = logrank_p, HR = hr, HR_lo = hr_lo, HR_hi = hr_hi)
}

# Per-facet-cell N / event counts for annotating km_facet() grids as a companion
# table (printed under the plot — ggsurvplot_facet is not a safe geom_text target).
# BASE R (no dplyr) so it is unit-testable in the bootstrap env. Reuses recode_event().
km_panel_counts <- function(df, facets,
                            time  = attend_cols$surv_time,
                            event = attend_cols$surv_event,
                            min_n = 10, min_events = 5) {
  keep <- c(time, event, facets)
  if (!all(keep %in% names(df))) return(data.frame())
  ev <- recode_event(df[[event]])
  ok <- !is.na(df[[time]]) & !is.na(ev)
  for (f in facets) ok <- ok & !is.na(df[[f]])
  if (!any(ok)) return(data.frame())
  g  <- lapply(facets, function(f) as.character(df[[f]][ok]))
  names(g) <- facets
  ev <- ev[ok]
  agg <- stats::aggregate(ev, by = g,
                          FUN = function(x) c(n = length(x), n_events = sum(x == 1)))
  out <- data.frame(agg[facets], n = agg$x[, "n"], n_events = agg$x[, "n_events"],
                    row.names = NULL, stringsAsFactors = FALSE)
  out$small_n <- out$n < min_n | out$n_events < min_events
  out
}

# =============================================================================
# --- 8. Copy-number / CNV configuration (report 08) -------------------------
# Reproduces the TCGA endometrial integrated classification (Kandoth et al.,
# Nature 2013; papers/nature12113.pdf). The CN-low vs CN-high split comes from
# clustering copy-number profiles; the high-aneuploidy cluster = "copy-number
# high (serous-like)". Two CNV inputs are supported, preferred in this order:
#
#   ASCETS  (Spurr et al., Genome Med 2021) — the BETTER path. Turns segmented
#           copy number (e.g. DRAGEN .seg) into a validated arm-level matrix +
#           a per-sample aneuploidy score, benchmarked against TCGA arm calls.
#           We cluster its `arm_weighted_average_segmeans` matrix directly.
#   ClassifyCNV  (variantalker `*.cnv.annotated.tsv`) — the FALLBACK. Per-sample
#           altered-segment calls (DUP/DEL + logRatio). We aggregate `logRatio`
#           to chromosome arms ourselves (segments_to_arm_matrix), treating
#           unlisted regions as copy-neutral (0).
# Everything is config-driven and knit-safe: absent inputs -> the report skips
# the CNV sections, mirroring have_cols() / requireNamespace().
# --- Copy-number clustering parameters (report 08) --------------------------
#
# These four live here, not at the top of the report that uses them, for the reason every
# other threshold does (gotcha #2): a method parameter referenced from prose in one file
# and defined in another resolves to nothing. 00-methods.Rmd discusses the cut; report 08
# applies it; both source this file, so both can name it.
#
# 1 - Pearson correlation distance with Ward.D2 linkage. TCGA Suppl. Methods S2 leaves the
# *copy-number* clustering metric unstated. The pipeline used 1 - Pearson on the grounds that
# S2 names that pair for its mRNA and METHYLATION clusterings — but it never names it for copy
# number, and on this feature space it is measurably the wrong choice.
#
# WHY EUCLIDEAN. Correlation distance is UNDEFINED for a flat profile, and a copy-number-quiet
# tumour is flat at every significant peak, so cluster_arm_matrix() must drop those samples
# rather than invent a distance for them. TCGA's own published cluster 1 IS the quiet cluster
# (mean |log2| 0.0043, against 0.216 for cluster 4), so correlation discards most of the one
# cluster it would most need to recover. Measured on TCGA UCEC 2013 (365 tumours, their hg19
# segments, their published 79 peaks), both burden-ordered and scored against CNA_CLUSTER_K4:
#
#   1 - Pearson   281/365 clustered    7 of 86 cluster-1 retained   ARI 0.222   41.1% diagonal
#   Euclidean     365/365 clustered   86 of 86 cluster-1 retained   ARI 0.423   64.7% diagonal
#
# Euclidean loses no tumours, recovers the quiet cluster completely, nearly doubles the ARI, and
# its contingency is clean on the diagonal (published 4 -> recomputed 4 at 79/93, no leakage
# into 1-3). A flat sample is not undefined under Euclidean; it simply sits near the origin,
# which is the correct geometry for "no alteration". Report 10 prints the comparison on every
# knit, so this is a claim the site re-checks rather than one it asserts.
#
# Still ONE distance/linkage combination for the whole pipeline — that rule is unchanged.
cnv_method <- "ward.D2"       # Ward linkage
cnv_dist   <- "euclidean"     # see above: correlation cannot represent a flat SCNA profile

# Cut height for the cascade / CN-high-burden label. TCGA found 4 SCNA clusters (Fig. 1a);
# ATTEND has no POLE and is smaller, but k = 4 keeps parity. Change freely.
k_cnv <- 4

# The chromosome heatmap is drawn ONCE with a cluster bar per k in this range, so k = 2..8
# can be compared in one figure (fig1a_heatmap_ksweep()).
k_sweep <- 2:8

attend_cnv <- list(
  build     = "hg38",                    # which arm-boundary table to use
  arm_table = "arm_boundaries_hg38.tsv", # bundled in data/
  # Per-assembly tables, selected by load_arm_boundaries(build=). ATTEND's own DRAGEN
  # .seg is hg38; the TCGA 2013 SNP6 segments are hg19 — same arm definitions,
  # different coordinates, so they must NOT share a boundary table.
  arm_tables = list(hg38 = "arm_boundaries_hg38.tsv",
                    hg19 = "arm_boundaries_hg19.tsv"),

  # ClassifyCNV-annotated segment files (variantalker output).
  classifycnv = list(
    dir       = "cnv_annotations",       # data/cnv_annotations/<s>/<s>.cnv.annotated.tsv
    chrom_col = "Chromosome",
    start_col = "Start",
    end_col   = "End",
    type_col  = "Type",
    value_col = "logRatio",              # primary clustering signal (log2 ratio)
    cnf_col   = "CNF",                    # 2^logRatio (linear fold)
    id_strip  = "-1TAD104|_tumor_only"   # normalise sequencing id -> bare barcode
  ),

  # ASCETS outputs (run DRAGEN .seg -> ASCETS, drop the *_*.txt files here).
  ascets = list(
    dir             = "ascets",          # data/ascets/
    armmeans_glob   = "*arm_weighted_average_segmeans*.txt",
    aneuploidy_glob = "*aneuploidy_scores*.txt",
    calls_glob      = "*arm_level_calls*.txt",
    id_strip        = "-1TAD104|_tumor_only"
  ),

  # DRAGEN genome-wide segmentation (.seg) — the input to the TCGA Fig-1a-style
  # GENOME-BINNED clustering (segments_to_bin_matrix()), and the same file that feeds
  # ASCETS. One .seg per sample under data/seg/; ID comes from the file name. Column
  # names are AUTO-DETECTED across DRAGEN variants (first match in each vector wins),
  # so unusual headers just need a name added here. `value_is_log2` = the seg-mean
  # column is a log2 ratio (DRAGEN default); set FALSE for a linear copy ratio / CN.
  seg = list(
    dir           = "seg",               # data/seg/  (one <sample>.seg per file; subfolders OK)
    chrom_col     = c("Chromosome", "chrom", "chr", "#chrom", "seqnames", "Chr"),
    start_col     = c("Start", "start", "loc.start", "Start_Position", "startpos"),
    end_col       = c("End", "end", "loc.end", "End_Position", "endpos"),
    value_col     = c("Segment_Mean", "seg.mean", "seg_mean", "SegmentMean",
                      "copy_ratio", "MeanLog2R", "log2", "Mean", "score"),
    num_col       = c("num.mark", "Num_Markers", "num_mark", "num.probes",
                      "Num_Probes", "N.BAF", "probes"),   # marker count (for GISTIC)
    value_is_log2 = TRUE,
    bin           = 1e6,                  # genome bin width (bp) for the clustering matrix
    id_strip      = "-1TAD104|_tumor_only",
    # GISTIC2 input: the per-sample .seg files are concatenated into ONE seg file here
    # by write_gistic_seg(); feed that to gistic2 (see code/run_gistic.sh). GISTIC runs
    # OUTSIDE R; its outputs land in data/gistic/ and report 06 and 35 pick them up.
    gistic_seg_out = "output/gistic_input/attend_all_segments.seg"
  ),

  # Report 07 oncoplot: how many of the MOST-recurrently-altered chromosome arms
  # (by cohort alteration frequency, from the ASCETS arm-level calls) to add as
  # annotation tracks alongside the mutation matrix. Keep small for a legible figure.
  oncoplot_arms = 5,

  # GISTIC2 focal peaks — the paper-standard FOCAL (gene-level) recurrent-CNV layer,
  # overlaid amp/del in the oncoplot grid (TCGA nature12113 style). Run GISTIC2 on the
  # DRAGEN .seg OUTSIDE R and drop its output in data/gistic/; report 06 then attaches
  # it to read.maf automatically. Absent -> the oncoplot is mutation+arm only (knit-safe).
  gistic = list(
    dir               = "gistic",                # data/gistic/
    all_lesions_glob  = "*all_lesions*.txt",     # all_lesions.conf_XX.txt
    amp_genes_glob    = "*amp_genes*.txt",       # amp_genes.conf_XX.txt
    del_genes_glob    = "*del_genes*.txt",       # del_genes.conf_XX.txt
    scores_glob       = "*scores.gistic",        # scores.gistic — REQUIRED by
                                                 # maftools::readGistic (read.maf forwards
                                                 # to it & hard-stops without this file)
    # THE faithful TCGA clustering input (Kandoth Suppl. Methods S2): "thresholded relative
    # copy number data in significantly reoccurring amp/del regions identified by GISTIC 2.0".
    # all_thresholded.by_genes.txt = genes x samples in {-2,-1,0,1,2}; we restrict it to
    # genes inside the significant peaks (amp_genes/del_genes) and cluster THAT.
    all_thresholded_glob = "*all_thresholded*.txt",

    # THE faithful TCGA *FIGURE* INPUT, which is a DIFFERENT MATRIX from the clustering one.
    # Fig. 1a is "SCNAs in each tumour (horizontal axis) plotted by chromosomal location
    # (vertical axis)" — the whole genome, on a continuous colour scale. The clustering input
    # above is neither: it is peak-restricted (most of the genome is absent) and discretised
    # to five values. Plotting the clustering matrix is what made report 09's heatmap not the
    # paper's heatmap. all_data_by_genes.txt is the same GISTIC run's genome-wide CONTINUOUS
    # copy number per gene, with the Cytoband column that places each row on the genome.
    # Cluster on the peaks, DISPLAY this — exactly as the paper does.
    all_data_glob = "*all_data_by_genes*.txt"
  ),

  # Recurrent chromosome-arm SCNAs (recurrent_arm_calls(), report 06): an arm must be
  # EVALUABLE (call != LOWCOV/NC/NA) in at least this fraction of samples to be eligible
  # for top_arms — keeps a poorly-covered arm from ranking on a handful of confident calls.
  # n_barplot: how many of the eligible, most-altered arms the report-5 barplot draws
  # (the axis is genome-ordered afterwards, so this caps the figure, not the eligibility
  # gate — the same `freq$eligible` gate applies before this slice, GC3/config not literal).
  arm_calls = list(min_evaluable_frac = 0.5, n_barplot = 20),

  # Genome-wide recurrent-SCNA pileup (report 06, TCGA nature12113 Fig-1 style):
  # base-pair-resolution fraction of samples gained/lost along the genome, from the
  # ClassifyCNV altered segments. `bin` = window width (bp). Direction comes from the
  # segment Type (DUP/DEL); min_abs_logratio can additionally require a magnitude.
  # high_abs_logratio: the low-level/high-level amplitude boundary (Beroukhim 2010; Zack
  # 2013; GenVisR::cnFreq separate "gain/loss" from "amplification/homozygous deletion"
  # the same way). |log2 ratio| >= 1 is roughly a full copy gained or lost relative to the
  # segment's baseline (log2(3/2) ~= 0.58 for a single-copy gain in a diploid tumour, so
  # >=1 requires more than that) — the conventional cutoff used to call amplitude tiers.
  pileup = list(bin = 1e6, min_abs_logratio = 0, high_abs_logratio = 1.0)
)

# Recurrent chromosome-arm SCNAs from an ASCETS arm-level calls table (ID + one column
# per arm; values AMP/DEL/NEUTRAL/NC, or +1/-1/0). PURE transform, no I/O — feed it
# load_ascets_calls(). Returns a list:
#   long     : ID, arm, call ("gain"/"loss"/"neutral"/NA) — tidy per-sample-per-arm
#   freq     : arm + cohort fractions (gain, loss, altered) + n_evaluable/n_lowcov +
#              eligible, sorted most-altered first. Denominator is EVALUABLE arms
#              (LOWCOV/NC/NA/"" excluded by norm_call()) — this matches ASCETS's own
#              aneuploidy score and load_wes_results.R:507's n_arms_evaluable ONLY on
#              LOWCOV. That line is `sum(!is.na(call) & call != "LOWCOV")`, which counts
#              NC as evaluable; here NC is a no-call exactly like LOWCOV (deliberate — see
#              norm_call() below), so the two denominators can differ when NC calls are
#              present. `eligible` is TRUE when an arm is evaluable in at least
#              `min_evaluable_frac` of samples — the SAME rule `top_arms` uses below, so a
#              caller (e.g. the report's barplot) can apply the identical gate to a
#              wider/differently-ranked slice of `freq` without recomputing it and risking
#              disagreement with `top_arms`.
#   top_arms : the `n_top` most-altered arm names among `freq$eligible` arms (the ones
#              worth annotating)
#   wide_top : ID + one "CNV_<arm>" column per top arm — the oncoplot annotation table
# NULL in -> NULL out (knit-safe).
recurrent_arm_calls <- function(calls, n_top = 5,
                                min_evaluable_frac = attend_cnv$arm_calls$min_evaluable_frac) {
  if (is.null(calls) || nrow(calls) == 0 || !"ID" %in% names(calls)) return(NULL)
  arm_cols <- setdiff(names(calls), "ID")
  if (length(arm_cols) == 0) return(NULL)
  # min_evaluable_frac denominator: DISTINCT samples in the input, not nrow(calls) — a
  # resequenced sample appearing twice (load_ascets_calls() applies the same id_strip
  # collapse used elsewhere in the pipeline, so this can happen) must not double-count,
  # consistent with how load_cnv_data()'s n_profiled is derived (distinct ids, not rows).
  n_id <- dplyr::n_distinct(calls$ID)
  norm_call <- function(x) {
    # x is already upper-cased + trimmed. "NA" the literal string and NA the real value
    # are the same no-call. LOWCOV/NC/""/NA all map to NA_character_ — they must NOT fall
    # into a "neutral" catch-all, or ASCETS's low-coverage flag silently becomes a real
    # call and inflates both the gain/loss numerators and the arm denominator.
    out <- dplyr::case_when(
      x %in% c("AMP", "GAIN", "1", "+1", "2")        ~ "gain",
      x %in% c("DEL", "LOSS", "-1", "-2")             ~ "loss",
      x %in% c("NEUTRAL", "NEUT", "0")                ~ "neutral",
      x %in% c("LOWCOV", "NC", "NA", "") | is.na(x)   ~ NA_character_,
      TRUE                                             ~ NA_character_)
    unrecognised <- unique(x[is.na(out) & !(x %in% c("LOWCOV", "NC", "NA", "") | is.na(x))])
    if (length(unrecognised) > 0)
      warning("norm_call(): unrecognised token(s), treated as no-call: ",
               paste(unrecognised, collapse = ", "))
    out
  }
  long <- calls |>
    pivot_longer(all_of(arm_cols), names_to = "arm", values_to = "call") |>
    mutate(call = norm_call(toupper(trimws(as.character(call)))))
  freq <- long |>
    group_by(arm) |>
    summarise(gain        = mean(call == "gain",    na.rm = TRUE),
              loss        = mean(call == "loss",    na.rm = TRUE),
              altered     = mean(call != "neutral", na.rm = TRUE),
              n_evaluable = sum(!is.na(call)),
              n_lowcov    = sum(is.na(call)),
              .groups = "drop") |>
    mutate(eligible = n_evaluable / n_id >= min_evaluable_frac) |>
    arrange(desc(altered))
  top_arms <- head(freq$arm[freq$eligible], n_top)
  wide_top <- long |>
    filter(arm %in% top_arms) |>
    mutate(arm = paste0("CNV_", arm)) |>
    pivot_wider(names_from = arm, values_from = call, values_fn = ~ .x[1])
  list(long = long, freq = freq, top_arms = top_arms, wide_top = wide_top)
}

# Order chromosome-arm names ALONG THE GENOME (1p, 1q, 2p, 2q, ..., 22p, 22q, Xp, Xq, Yp,
# Yq) rather than alphabetically — field convention for recurrent-CNA figures (Beroukhim
# et al. 2010; Zack et al. 2013; GenVisR::cnFreq) is to display recurrence along the
# genome, and a plain string sort gets this wrong ("10p" sorts before "2p", and even before
# "1q"). This is a PURE reordering: every input element appears exactly once in the output,
# nothing is invented or dropped, and it does NOT touch recurrent_arm_calls()'s own
# most-altered-first sort (freq/top_arms) — the caller decides which order to display in
# (GC4). Tolerates a "chr" prefix and mixed case ("chr1p", "CHR1P", "1P" all parse the
# same). Names that don't match `<chr>?<1-22|X|Y><p|q>` sort LAST, in their original
# relative order (a stable sort on an unparseable key), rather than erroring.
order_arms_genomic <- function(arms) {
  chrom_rank <- c(as.character(1:22), "X", "Y")
  parsed <- stringr::str_match(toupper(trimws(arms)), "^(?:CHR)?([0-9]+|X|Y)([PQ])$")
  chrom_key <- match(parsed[, 2], chrom_rank)      # NA when the chromosome part is unparseable
  arm_key   <- match(parsed[, 3], c("P", "Q"))     # NA when the p/q part is unparseable
  key <- chrom_key * 2 + arm_key                   # NA propagates -> unparseable names sort last
  arms[order(is.na(key), key, seq_along(arms))]
}

# Genome-wide recurrent-SCNA PILEUP (report 06). PURE transform, no I/O — feed it
# load_cnv_data() (ClassifyCNV altered segments) and load_arm_boundaries(). Tiles the
# genome into `bin`-wide windows and, per window, returns the fraction of samples with
# a gain / loss segment overlapping it — the base-pair companion to the arm-level
# recurrence barplot (TCGA nature12113 Fig. 1). Direction is taken from the segment
# `Type` (DUP/DEL) with the `logRatio` sign as a fallback; `min_abs` optionally drops
# small-magnitude segments. Chromosome lengths + cumulative genome offsets come from
# the arm-boundary table (q-arm end per chromosome). Returns a tibble
#   chrom, bin_start, gpos (cumulative genome coordinate), gain, loss,
#   gain_low, gain_high, loss_low, loss_high
# `gain`/`loss` keep their original meaning — fraction with ANY gain/loss segment
# overlapping the bin (GC4, unchanged). `gain_low`/`gain_high` (and the loss pair) split
# that same population into low- vs high-amplitude tiers by `|logRatio| >= high_abs`
# (field convention: low-level gain/loss vs high-level amplification/homozygous deletion —
# Beroukhim 2010; Zack 2013; GenVisR::cnFreq). Tiering is per (bin, direction, sample): a
# sample with both a low- and a high-magnitude segment overlapping the same bin in the same
# direction is counted ONCE, in the high tier only, so `gain_low + gain_high == gain`
# (and likewise for loss) always. `NA` logRatio counts as low tier.
# with attr "chrom_offsets" (for axis labelling) and "n_samples". NULL in -> NULL out.
segment_pileup <- function(cnv_long, arms, bin = attend_cnv$pileup$bin,
                           min_abs = attend_cnv$pileup$min_abs_logratio,
                           chrom_col = "Chromosome", start_col = "Start",
                           end_col = "End", type_col = "Type",
                           value_col = "logRatio", id_col = "ID",
                           n_samples = NULL,
                           high_abs = attend_cnv$pileup$high_abs_logratio) {
  if (is.null(cnv_long) || nrow(cnv_long) == 0 || is.null(arms)) return(NULL)

  # Denominator = samples PROFILED, not samples with a surviving altered segment. Three
  # sources, tried in this order, recorded in attr(wide, "n_samples_source") so the report
  # can state its denominator honestly:
  #   "argument" — an explicit n_samples always wins (e.g. the full cohort size, covering
  #                samples whose annotation file failed to parse entirely).
  #   "profiled" — attr(cnv_long, "n_profiled"), set by load_cnv_data() from the count of
  #                DISTINCT derived sample ids among the per-sample files it read (not raw
  #                file count — two files can de-duplicate to one sample id, e.g. a resequence
  #                under a different id_strip suffix). This is the ONLY source that can see a
  #                sample that contributed ZERO rows to cnv_long (a copy-number-quiet tumour whose
  #                ClassifyCNV file has no calls) — no inspection of cnv_long's rows can
  #                recover that count. Used only when present, non-NA, numeric, and at
  #                least the observed distinct-id count below (a floor guard: a stale or
  #                wrong attribute can never SHRINK the denominator below the samples
  #                actually in hand).
  #   "observed" — dplyr::n_distinct() on the UNFILTERED cnv_long, before min_abs /
  #                direction drop rows (Task 1) — a sample whose genome is
  #                quiet-but-present (or whose segments are all below min_abs) still
  #                counts.
  n_observed <- dplyr::n_distinct(cnv_long[[id_col]])
  if (!is.null(n_samples)) {
    n_samples_source <- "argument"
  } else {
    profiled <- attr(cnv_long, "n_profiled")
    if (!is.null(profiled) && length(profiled) == 1 && !is.na(profiled) &&
        is.numeric(profiled) && profiled >= n_observed) {
      n_samples        <- profiled
      n_samples_source <- "profiled"
    } else {
      n_samples        <- n_observed
      n_samples_source <- "observed"
    }
  }

  segs <- cnv_long |>
    transmute(
      ID    = as.character(.data[[id_col]]),
      chrom = as.character(.data[[chrom_col]]),
      s     = suppressWarnings(as.numeric(.data[[start_col]])),
      e     = suppressWarnings(as.numeric(.data[[end_col]])),
      type  = if (type_col  %in% names(cnv_long)) toupper(trimws(as.character(.data[[type_col]]))) else NA_character_,
      v     = if (value_col %in% names(cnv_long)) suppressWarnings(as.numeric(.data[[value_col]]))    else NA_real_) |>
    filter(!is.na(ID), !is.na(chrom), !is.na(s), !is.na(e), e > s)
  if (!is.null(min_abs) && min_abs > 0)
    segs <- segs |> filter(is.na(v) | abs(v) >= min_abs)
  if (nrow(segs) == 0) return(NULL)

  segs <- segs |>
    mutate(direction = case_when(
      type %in% c("DUP", "GAIN", "AMP") ~ "gain",
      type %in% c("DEL", "LOSS")        ~ "loss",
      !is.na(v) & v > 0                 ~ "gain",
      !is.na(v) & v < 0                 ~ "loss",
      TRUE                              ~ NA_character_),
      # amplitude tier per SEGMENT; NA logRatio counts as low (ambiguity resolution #5).
      # Tier is later collapsed per (bin, direction, sample) so a sample never double-counts.
      tier = ifelse(!is.na(v) & abs(v) >= high_abs, "high", "low")) |>
    filter(!is.na(direction))
  if (nrow(segs) == 0) return(NULL)

  # chromosome lengths (q-arm end) + cumulative offsets, genome order 1..22, X, Y
  chrom_len <- arms |>
    transmute(chrom = as.character(chrom), end = as.numeric(end)) |>
    group_by(chrom) |> summarise(len = max(end), .groups = "drop") |>
    mutate(ord = match(str_remove(chrom, "^chr"), c(as.character(1:22), "X", "Y"))) |>
    filter(!is.na(ord)) |>
    arrange(ord) |>
    mutate(offset = cumsum(dplyr::lag(len, default = 0)),
           mid    = offset + len / 2)
  segs <- segs |> filter(chrom %in% chrom_len$chrom)
  if (nrow(segs) == 0) return(NULL)

  # tile each chromosome into `bin`-wide windows
  bins <- chrom_len |>
    transmute(chrom, len, offset,
              bin_start = purrr::map(len, ~ seq(0, .x - 1, by = bin))) |>
    tidyr::unnest(bin_start) |>
    mutate(bin_end = pmin(bin_start + bin, len),
           gpos    = offset + bin_start + bin / 2)

  # segments overlapping each bin — the shared join both "gain"/"loss" (GC4, any-tier)
  # and the amplitude-tier columns are built from, so the two can never disagree about
  # which segments overlap which bin.
  bin_segs <- bins |>
    select(chrom, bin_start, bin_end) |>
    inner_join(segs, by = "chrom", relationship = "many-to-many") |>
    filter(s < bin_end, e > bin_start)

  # per bin, fraction of DISTINCT samples with a gain / loss segment overlapping it
  # (unchanged meaning from before amplitude tiers were added — GC4).
  overlap <- bin_segs |>
    distinct(chrom, bin_start, direction, ID) |>
    count(chrom, bin_start, direction, name = "n_alt") |>
    mutate(frac = n_alt / n_samples)

  # amplitude tier per (bin, direction, sample): a sample with both a low- and a
  # high-magnitude segment overlapping the same bin+direction counts ONCE, in the high
  # tier only (ambiguity resolution #4) — so tier_pick partitions the same distinct-ID
  # population `overlap` counts, and gain_low + gain_high == gain (and likewise loss).
  overlap_tier <- bin_segs |>
    group_by(chrom, bin_start, direction, ID) |>
    summarise(tier_pick = if (any(tier == "high")) "high" else "low", .groups = "drop") |>
    count(chrom, bin_start, direction, tier_pick, name = "n_alt") |>
    mutate(frac = n_alt / n_samples)

  tier_col <- function(dir, tier_val)
    overlap_tier |>
      filter(direction == dir, tier_pick == tier_val) |>
      select(chrom, bin_start, frac)

  wide <- bins |>
    select(chrom, bin_start, gpos) |>
    left_join(overlap |> filter(direction == "gain") |> select(chrom, bin_start, gain = frac),
              by = c("chrom", "bin_start")) |>
    left_join(overlap |> filter(direction == "loss") |> select(chrom, bin_start, loss = frac),
              by = c("chrom", "bin_start")) |>
    left_join(tier_col("gain", "low")  |> rename(gain_low  = frac), by = c("chrom", "bin_start")) |>
    left_join(tier_col("gain", "high") |> rename(gain_high = frac), by = c("chrom", "bin_start")) |>
    left_join(tier_col("loss", "low")  |> rename(loss_low  = frac), by = c("chrom", "bin_start")) |>
    left_join(tier_col("loss", "high") |> rename(loss_high = frac), by = c("chrom", "bin_start")) |>
    mutate(across(c(gain, loss, gain_low, gain_high, loss_low, loss_high),
                  ~ tidyr::replace_na(.x, 0)))

  attr(wide, "chrom_offsets")    <- chrom_len
  attr(wide, "n_samples")        <- n_samples
  attr(wide, "n_samples_source") <- n_samples_source
  wide
}

# --- 9. POLE-ultramutated configuration (NO LONGER in the cascade) ----------
# ATTEND has no POLE-ultramutated cases, so the POLE branch has been REMOVED from
# the integrated cascade (add_tcga_class starts at MSI). pole_ultramutated() and
# this config are kept, unwired, in case a POLE branch is wanted again — call
# pole_ultramutated() and reinstate the branch in add_tcga_class() to re-enable.
# Calling is by exonuclease-domain hotspot (needs a protein-change column in the
# MAF) OR, as a backstop, an ultra-high TMB cutoff.
attend_pole <- list(
  gene        = "POLE",
  hotspots    = c("P286R", "V411L", "S297F", "A456P", "S459F", "F367S", "L424V"),
  protein_col = attend_maf$protein_col,   # one candidate list, defined in attend_maf
  tmb_cutoff  = 100                      # mut/Mb; ultra-high backstop
)

# --- 10. Integrated-class labels (report 08) --------------------------------
# POLE is intentionally omitted: ATTEND has no POLE-ultramutated cases and the
# cascade starts at MSI. (pole_ultramutated()/attend_pole are kept below, unwired,
# in case a POLE branch is wanted again.)
attend_tcga_levels <- list(
  # Top cascade tier. Criterion is MMR-deficiency by IHC ALONE (MMRd vs MMRp) — see
  # is_mmrd()/add_tcga_class(). MSI-high (NGS) is a correlate, NOT a criterion. The list
  # key stays `msi` for back-compat with existing references; the displayed label is "MMRd".
  msi     = "MMRd",
  cn_low  = "Copy-number low (endometrioid)",
  cn_high = "Copy-number high (serous-like)"
)

# How the copy-number-high (serous-like) call is made among the POLE-wt / MMR-proficient
# samples. TCGA's LITERAL definition is SCNA cluster-4 membership (Kandoth et al. Fig. 2b:
# "serous-like ... copy-number cluster 4") — TP53 is a *property* of that cluster (~90%
# mutant), NOT a criterion. So the FAITHFUL reproduction is "cluster_only"; TP53 is then
# validated separately (report 08 "TP53 enrichment" QC), not required. Set here:
#   "cluster_only" serous = high-aneuploidy CNV cluster                       (DEFAULT; faithful TCGA repro)
#   "and"          serous = high-aneuploidy CNV cluster AND TP53-pathogenic   (higher-precision variant)
#   "or"           serous = high-aneuploidy CNV cluster OR  TP53-pathogenic   (workflow §1.1 union rule)
#   "tp53_only"    serous = TP53-pathogenic                                   (portable, no clustering)
# If TP53 status is absent, "and"/"tp53_only" degrade to "cluster_only" (knit-safe).
attend_tcga_serous_rule <- "cluster_only"

# --- 11. CNV helpers --------------------------------------------------------

# Load the bundled chromosome-arm boundary table (chrom, arm, start, end).
# Returns NULL (with a message) if missing, so the report degrades gracefully.
# `build` selects the assembly: NULL (default) keeps the pipeline build (attend_cnv$build,
# hg38 — ATTEND's DRAGEN .seg), "hg19" the bundled hg19 table. The hg19 table exists for
# note: the TCGA 2013 SNP6 segments are hg19, and binning hg19 coordinates against
# hg38 chromosome lengths would shift every genome bin and misplace the centromere.
load_arm_boundaries <- function(cnv = attend_cnv, build = NULL) {
  tbl <- if (is.null(build)) cnv$arm_table
         else if (!is.null(cnv$arm_tables[[build]])) cnv$arm_tables[[build]]
         else { message("No arm table for build '", build, "' — using ", cnv$arm_table); cnv$arm_table }
  p <- here::here("data", tbl)
  if (!file.exists(p)) { message("Arm-boundary table not found: ", p); return(NULL) }
  readr::read_tsv(p, show_col_types = FALSE)
}

# Aggregate long CNV segments to a sample x arm matrix of length-weighted mean
# copy-number. Each segment is clipped to the arms it overlaps; the per-arm value
# is sum(overlap_len * value) / arm_length, so regions NOT listed in the file
# contribute 0 (copy-neutral) — making the value a magnitude-weighted
# fraction-of-arm-altered (a broad-SCNA-burden measure). Returns a tibble:
# ID + one numeric column per arm (missing arms filled 0). Empty input -> NULL.
segments_to_arm_matrix <- function(cnv_long, arms,
                                   value_col = attend_cnv$classifycnv$value_col,
                                   chrom_col = "Chromosome",
                                   start_col = "Start",
                                   end_col   = "End",
                                   id_col    = "ID") {
  if (is.null(cnv_long) || nrow(cnv_long) == 0 || is.null(arms)) return(NULL)

  segs <- cnv_long |>
    transmute(ID    = as.character(.data[[id_col]]),
              chrom = as.character(.data[[chrom_col]]),
              s     = suppressWarnings(as.numeric(.data[[start_col]])),
              e     = suppressWarnings(as.numeric(.data[[end_col]])),
              v     = suppressWarnings(as.numeric(.data[[value_col]]))) |>
    filter(!is.na(ID), !is.na(chrom), !is.na(s), !is.na(e), !is.na(v))
  if (nrow(segs) == 0) return(NULL)

  arms2 <- arms |>
    transmute(chrom   = as.character(chrom),
              arm     = as.character(arm),
              a_s     = as.numeric(start),
              a_e     = as.numeric(end),
              arm_len = as.numeric(end) - as.numeric(start))
  arm_levels <- arms2$arm

  # Align segment chrom naming to the arm table's convention (chr-prefixed or bare),
  # exactly as segments_to_bin_matrix() does — otherwise the join silently yields
  # nothing when the .seg uses "chr1" and the table "1" (or vice versa) and the arm
  # matrix comes back empty.
  bt_has_chr <- all(str_detect(arms2$chrom, "^chr"))
  segs <- segs |> mutate(chrom = if (bt_has_chr) paste0("chr", str_remove(chrom, "^chr"))
                                 else str_remove(chrom, "^chr"))

  segs |>
    inner_join(arms2, by = "chrom", relationship = "many-to-many") |>
    mutate(ov = pmax(0, pmin(e, a_e) - pmax(s, a_s))) |>
    filter(ov > 0) |>
    group_by(ID, arm, arm_len) |>
    summarise(weighted = sum(ov * v), .groups = "drop") |>
    transmute(ID, arm, arm_value = weighted / arm_len) |>
    tidyr::complete(ID, arm = arm_levels, fill = list(arm_value = 0)) |>
    pivot_wider(names_from = arm, values_from = arm_value, values_fill = 0) |>
    select(ID, any_of(arm_levels))   # stable arm column order
}

# Genome-BINNED per-sample copy-number matrix — the faithful TCGA Fig-1a clustering
# feature (Kandoth et al. clustered SNP6 segmented CN across genome position, not arms).
# Tiles each chromosome into `bin`-wide windows and, per sample x bin, takes the
# length-weighted mean seg-mean (log2 ratio) of the DRAGEN .seg segments overlapping it;
# genome not covered by a segment = 0 (copy-neutral). `seg` is the long table from
# load_seg_data() (ID, Chromosome, Start, End, Segment_Mean); `arms` supplies chromosome
# lengths + genome order (load_arm_boundaries()). Returns ID + one numeric column per bin
# ("chr1:0", "chr1:1000000", ...), columns in genome order (so a heatmap with row-
# clustering OFF shows the chromosomal axis). NULL on empty input. Cluster the result
# with cluster_arm_matrix() exactly as for the arm matrix (it is generic over ID+numeric).
segments_to_bin_matrix <- function(seg, arms, bin = attend_cnv$seg$bin,
                                   chrom_col = "Chromosome", start_col = "Start",
                                   end_col = "End", value_col = "Segment_Mean",
                                   id_col = "ID") {
  if (is.null(seg) || nrow(seg) == 0 || is.null(arms)) return(NULL)

  chrom_len <- arms |>
    transmute(chrom = as.character(chrom), end = as.numeric(end)) |>
    group_by(chrom) |> summarise(len = max(end), .groups = "drop") |>
    mutate(ord = match(str_remove(chrom, "^chr"), c(as.character(1:22), "X", "Y"))) |>
    filter(!is.na(ord)) |> arrange(ord)

  bins <- chrom_len |>
    transmute(chrom, ord, bin_start = purrr::map(len, ~ seq(0, .x - 1, by = bin))) |>
    tidyr::unnest(bin_start) |>
    mutate(bin_end = bin_start + bin,
           binid   = paste0(chrom, ":", format(bin_start, scientific = FALSE, trim = TRUE))) |>
    arrange(ord, bin_start)
  bin_levels <- bins$binid                        # genome order for heatmap rows

  segs <- seg |>
    transmute(ID    = as.character(.data[[id_col]]),
              chrom = as.character(.data[[chrom_col]]),
              s     = suppressWarnings(as.numeric(.data[[start_col]])),
              e     = suppressWarnings(as.numeric(.data[[end_col]])),
              v     = suppressWarnings(as.numeric(.data[[value_col]]))) |>
    filter(!is.na(ID), !is.na(chrom), !is.na(s), !is.na(e), !is.na(v), e > s)
  # normalise chrom to match the boundary table's "chr" prefixing
  bt_has_chr <- all(str_detect(chrom_len$chrom, "^chr"))
  segs <- segs |> mutate(chrom = if (bt_has_chr) paste0("chr", str_remove(chrom, "^chr"))
                                 else str_remove(chrom, "^chr"))
  if (nrow(segs) == 0) return(NULL)

  bins |>
    select(chrom, bin_start, bin_end, binid) |>
    inner_join(segs, by = "chrom", relationship = "many-to-many") |>
    mutate(ov = pmax(0, pmin(e, bin_end) - pmax(s, bin_start))) |>
    filter(ov > 0) |>
    group_by(ID, binid) |>
    summarise(v = sum(ov * v) / bin, .groups = "drop") |>        # length-weighted log2
    tidyr::complete(ID, binid = bin_levels, fill = list(v = 0)) |>
    pivot_wider(names_from = binid, values_from = v, values_fill = 0) |>
    select(ID, any_of(bin_levels))                               # stable genome order
}

# Hierarchical clustering of a sample x feature matrix (arm or genome-bin) —
# matching TCGA's "unsupervised hierarchical clustering". `arm_tbl` is ID + feature
# columns (segments_to_arm_matrix() / segments_to_bin_matrix() / ASCETS arm-means /
# a GISTIC gene matrix). `dist_method` = the distance ("euclidean" default; set
# "correlation" for 1 - Pearson, which some TCGA CN clusterings use — see the paper's
# Supplementary Methods to match exactly) and `method` = the linkage (Ward.D2 default).
# Returns the numeric matrix (rownames = ID), the hclust, the cut vector, and a
# tidy ID->cnv_cluster assignment. NULL if too few samples for k clusters.
cluster_arm_matrix <- function(arm_tbl, k = 4, id_col = "ID",
                               method = "ward.D2", dist_method = "euclidean") {
  if (is.null(arm_tbl) || nrow(arm_tbl) < k) {
    message("cluster_arm_matrix(): need >= ", k, " samples, have ",
            if (is.null(arm_tbl)) 0 else nrow(arm_tbl), " — skipping."); return(NULL)
  }
  m <- as.data.frame(arm_tbl)
  rownames(m) <- m[[id_col]]; m[[id_col]] <- NULL
  m <- as.matrix(m)
  m <- m[, apply(m, 2, function(col) !all(is.na(col))), drop = FALSE]  # drop all-NA features
  m[is.na(m)] <- 0
  m <- m[stats::complete.cases(m), , drop = FALSE]
  if (nrow(m) < k) { message("cluster_arm_matrix(): < ", k, " complete rows."); return(NULL) }

  # 1 - Pearson correlation distance (across features) when requested, else a base dist metric.
  if (identical(dist_method, "correlation")) {
    # Correlation is UNDEFINED for a sample with zero variance across the features — and
    # such samples are common on thresholded GISTIC calls, where a near-diploid tumour can
    # be flat 0 at every peak gene. cor() returns NA for those rows and hclust() then dies
    # with "NA/NaN/Inf in foreign function call". Drop them explicitly and SAY SO: they
    # cannot be placed by a correlation metric, and inventing a distance for them would
    # silently invent structure. (Euclidean has no such problem — a flat sample is simply
    # near the origin — so the two dist_methods can legitimately cluster different sample
    # sets; downstream tables report n per clustering for exactly this reason.)
    sdv  <- apply(m, 1, stats::sd, na.rm = TRUE)
    flat <- !is.finite(sdv) | sdv == 0
    if (any(flat)) {
      message("cluster_arm_matrix(): dropping ", sum(flat), " zero-variance sample(s) — ",
              "correlation distance is undefined for a flat profile (", dist_method, ").")
      m <- m[!flat, , drop = FALSE]
      if (nrow(m) < k) {
        message("cluster_arm_matrix(): < ", k, " rows left after dropping flat profiles."); return(NULL)
      }
    }
    cm <- stats::cor(t(m), use = "pairwise.complete.obs")
    if (anyNA(cm)) {                       # backstop: treat any residual NA as uncorrelated
      message("cluster_arm_matrix(): ", sum(is.na(cm)), " undefined correlation(s) set to 0.")
      cm[is.na(cm)] <- 0
    }
    d <- stats::as.dist(1 - cm)
  } else {
    d <- stats::dist(m, method = dist_method)
  }
  hc <- stats::hclust(d, method = method)
  cl <- stats::cutree(hc, k = k)
  list(matrix     = m,
       hclust     = hc,
       cluster    = cl,
       assignment = tibble::tibble(ID = names(cl), cnv_cluster = unname(cl)))
}

# Adjusted Rand index between two clusterings of the same samples (Hubert & Arabie 1985).
# Chance-corrected agreement: 1 = identical partitions, ~0 = agreement no better than
# random, negative = worse than random. Cluster LABELS are ignored (only the partition
# matters), which is what makes it the right score for "did our recomputed clustering
# recover TCGA's published one" — their cluster 4 need not be our cluster 4.
# `a` and `b` are two label vectors of equal length; pairs where either is NA are dropped.
# Base R only (no mclust dependency). Returns NA_real_ for < 2 usable samples.
adjusted_rand_index <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  a <- as.character(a)[ok]; b <- as.character(b)[ok]
  n <- length(a)
  if (n < 2) return(NA_real_)
  tab   <- table(a, b)
  ch2   <- function(x) sum(x * (x - 1) / 2)             # choose(x, 2), vectorised
  idx   <- ch2(as.vector(tab))
  ra    <- ch2(rowSums(tab)); rb <- ch2(colSums(tab))
  total <- n * (n - 1) / 2
  expct <- ra * rb / total
  maxi  <- (ra + rb) / 2
  if (isTRUE(all.equal(maxi, expct))) return(NA_real_)  # degenerate (one cluster each)
  (idx - expct) / (maxi - expct)
}

# Per-cluster INTRINSIC SCNA burden, straight from the arm matrix: for each sample the
# mean |arm copy number| (magnitude of gains+losses across arms), averaged per cluster.
# This uses ONLY the clustering features, so it is INDEPENDENT of any external aneuploidy
# score — which is the whole point: label the copy-number-high cluster with THIS, then
# the aneuploidy score is free to serve as a non-circular enrichment check. `cl` is the
# cluster_arm_matrix() output (needs $matrix + $cluster). Returns a tibble
#   cnv_cluster, mean_burden, n  — sorted highest-burden first. NULL if unavailable.
#
# `direction`: NULL (default) reproduces the historical behaviour — mean |value|, which
# treats a deep deletion and a high amplification as equally "serous-like".
#
# ⚠️ WHY A DIRECTION ARGUMENT EXISTS. mean |value| answers "which cluster has the MOST
# alteration", which is NOT the same question as "which cluster is the serous-like one", and
# on THRESHOLDED features the two come apart. Measured on TCGA's published peaks:
#   cluster   n    %published-4   mean|value|   amp-gains   del-losses   directional
#      4     23        13%          0.829         0.310       0.569         0.879
#      3     79       100%          0.764         0.625       0.457         1.082
# The abs() rule picks cluster 4 — a small DELETION-dominated group — over cluster 3, which
# is 79/79 the published serous-like cluster. Supplying `direction` scores only CONCORDANT
# events (gains at amplification peaks, losses at deletion peaks) and recovers cluster 3.
# This is the same rule attend_scna already states for the report-15 panel score: "Panel loci
# are DIRECTIONAL: a deletion at MYC must not score as serous-like."
#
# `direction` is a NAMED character vector over the matrix's columns with values "amp"/"del"
# (anything else, or a missing name, falls back to |value| for that feature). Build it with
# gistic_feature_direction() for a GISTIC matrix, or from the published-peak table's
# `direction` column. Continuous feature spaces (genome bins, arm means) have no intrinsic
# per-feature direction — leave it NULL there.
cluster_scna_burden <- function(cl, direction = NULL) {
  if (is.null(cl) || is.null(cl$matrix) || is.null(cl$cluster)) return(NULL)
  ids <- rownames(cl$matrix)
  if (is.null(ids) || length(ids) == 0) return(NULL)
  M <- cl$matrix

  if (is.null(direction)) {
    burden <- rowMeans(abs(M), na.rm = TRUE)
  } else {
    d   <- direction[colnames(M)]                       # NA for unmapped features
    amp <- which(!is.na(d) & d == "amp")
    del <- which(!is.na(d) & d == "del")
    oth <- setdiff(seq_len(ncol(M)), c(amp, del))
    if (length(amp) + length(del) == 0) {
      message("cluster_scna_burden(): no feature matched the direction map — using |value|.")
      burden <- rowMeans(abs(M), na.rm = TRUE)
    } else {
      if (length(oth))
        message("cluster_scna_burden(): ", length(oth), " of ", ncol(M),
                " features have no direction — scored as |value|.")
      # concordant magnitude only: gains at amp loci, losses at del loci, |x| elsewhere
      contrib <- cbind(pmax(M[, amp, drop = FALSE],  0),
                       pmax(-M[, del, drop = FALSE], 0),
                       abs(M[, oth, drop = FALSE]))
      burden <- rowMeans(contrib, na.rm = TRUE)
    }
  }

  tibble::tibble(ID          = ids,
                 cnv_cluster = unname(cl$cluster[ids]),
                 scna_burden = burden[ids]) |>
    group_by(cnv_cluster) |>
    summarise(mean_burden = mean(scna_burden, na.rm = TRUE),
              n           = dplyr::n(), .groups = "drop") |>
    arrange(desc(mean_burden))
}

# Build the gene -> "amp"/"del" map for a GISTIC feature matrix, from that run's own
# amp_genes/del_genes peak lists. Genes appearing in BOTH lists are dropped (ambiguous), so
# they fall back to |value| rather than being assigned an arbitrary direction. Returns NULL
# when the peak lists are unavailable — callers then get the historical abs() behaviour.
gistic_feature_direction <- function(cnv = attend_cnv) {
  f <- tryCatch(find_gistic_files(cnv), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  a <- tryCatch(.gistic_peak_genes(f$amp_genes), error = function(e) character(0))
  d <- tryCatch(.gistic_peak_genes(f$del_genes), error = function(e) character(0))
  if (!length(a) && !length(d)) return(NULL)
  both <- intersect(a, d)
  a <- setdiff(a, both); d <- setdiff(d, both)
  if (length(both))
    message("gistic_feature_direction(): ", length(both),
            " gene(s) in both amp and del peaks — left undirected.")
  stats::setNames(c(rep("amp", length(a)), rep("del", length(d))), c(a, d))
}

# The copy-number-high / serous-like cluster = the cluster with the greatest intrinsic
# SCNA burden (cluster_scna_burden()). PREFERRED over cnv_high_cluster() because it does
# not consult the aneuploidy score, so the aneuploidy-by-cluster enrichment plot is a
# genuine independent check rather than a tautology. Returns the cluster id, or NA.
cnv_high_cluster_burden <- function(cl, direction = NULL) {
  b <- cluster_scna_burden(cl, direction = direction)
  if (is.null(b) || nrow(b) == 0) return(NA_integer_)
  b$cnv_cluster[[1]]
}

# Given an ID->cnv_cluster assignment and a per-ID aneuploidy value, return the
# cluster id with the highest mean aneuploidy = the copy-number-high / serous-
# like cluster. `aneu_tbl` must have columns ID + <aneu_col>. NOTE: this ranks by the
# SAME aneuploidy score you would then plot per cluster — so if that plot is your CN-high
# evidence, use cnv_high_cluster_burden() instead (independent signal) to avoid circularity.
cnv_high_cluster <- function(assignment, aneu_tbl, aneu_col = "aneuploidy_value") {
  if (is.null(assignment) || is.null(aneu_tbl)) return(NA_integer_)
  d <- assignment |> inner_join(aneu_tbl, by = "ID") |>
    filter(!is.na(.data[[aneu_col]]))
  if (nrow(d) == 0) return(NA_integer_)
  d |> group_by(cnv_cluster) |>
    summarise(m = mean(.data[[aneu_col]], na.rm = TRUE), .groups = "drop") |>
    slice_max(m, n = 1, with_ties = FALSE) |> pull(cnv_cluster)
}

# --- 12. POLE detection + integrated cascade --------------------------------

# Per-patient POLE-ultramutated flag. Prefers exonuclease-domain hotspots when a
# protein-change column is present in `maf_long`; otherwise falls back to an
# ultra-high TMB cutoff via `tmb_tbl` (ID + TMB column). Returns pid + logical.
# Expected to be all FALSE in ATTEND (no POLE-ultramutated cases).
pole_ultramutated <- function(maf_long = NULL, crosswalk = NULL, tmb_tbl = NULL,
                              tmb_col  = "TMB_SCORE",
                              pole_cfg = attend_pole) {
  out <- if (!is.null(crosswalk)) crosswalk |> distinct(pid) else tibble(pid = character())
  out$POLE_ultramutated <- FALSE
  if (nrow(out) == 0) return(out)

  prot_col <- if (!is.null(maf_long)) intersect(pole_cfg$protein_col, names(maf_long)) else character()
  if (length(prot_col) && !is.null(crosswalk)) {
    pat <- paste(pole_cfg$hotspots, collapse = "|")
    hit <- maf_long |>
      transmute(barcode = as.character(.data[[attend_maf$sample_col]]),
                gene    = as.character(.data[[attend_maf$gene_col]]),
                prot    = as.character(.data[[prot_col[1]]])) |>
      filter(gene == pole_cfg$gene, str_detect(replace_na(prot, ""), pat)) |>
      inner_join(crosswalk, by = "barcode") |>
      distinct(pid) |> mutate(.hot = TRUE)
    out <- out |> left_join(hit, by = "pid") |>
      mutate(POLE_ultramutated = POLE_ultramutated | replace_na(.hot, FALSE)) |>
      select(-.hot)
    return(out)
  }

  if (!is.null(tmb_tbl) && !is.null(crosswalk) && tmb_col %in% names(tmb_tbl)) {
    hit <- tmb_tbl |>
      transmute(barcode = as.character(ID), tmb = suppressWarnings(as.numeric(.data[[tmb_col]]))) |>
      inner_join(crosswalk, by = "barcode") |>
      group_by(pid) |> summarise(tmb = max(tmb, na.rm = TRUE), .groups = "drop") |>
      filter(is.finite(tmb), tmb >= pole_cfg$tmb_cutoff) |>
      distinct(pid) |> mutate(.hot = TRUE)
    out <- out |> left_join(hit, by = "pid") |>
      mutate(POLE_ultramutated = POLE_ultramutated | replace_na(.hot, FALSE)) |>
      select(-.hot)
  } else {
    message("pole_ultramutated(): no protein-change column and no TMB table — ",
            "POLE branch uncallable, all FALSE (expected in ATTEND).")
  }
  out
}

# Apply the TCGA integrated-classification PRIORITY CASCADE row-wise, POLE OMITTED:
#   MMRd (MSI-hypermutated)  >  copy-number-high (serous-like)  >  copy-number-low.
# Inputs are columns already on `df`: `msi_class_col` (== "MSI-high" from
# add_molecular_classes()), `mmr_class_col` (IHC MMR status, "Deficient"/"Intact"),
# `cnv_high_col` (logical, in the aneuploidy-enriched CNV cluster), and `tp53_col`
# (logical TP53_pathogenic).
#
# Top of the cascade — MMRd — is called when EITHER the NGS MSI caller is MSI-high OR
# MMR-protein loss is seen by IHC (workflow §1.1 step 2). MMRd overwrites the CN split.
#
# Serous-like (CN-high) among the MMR-proficient samples is set by `serous_rule`
# (attend_tcga_serous_rule): "cluster_only" (high-aneuploidy cluster membership alone =
# faithful TCGA repro; DEFAULT — TP53 is a *correlate* of the cluster, ~90% mutant,
# validated by report 08's TP53-enrichment QC, NOT a criterion), "and" (also require
# TP53-mut — higher-precision variant), "or" (the spec's union), or "tp53_only". Under the
# default "cluster_only", a high-cluster TP53-wild-type tumour is still called serous-like;
# only the opt-in "and" rule would demote it to CN-low. If TP53 status is absent,
# "and"/"tp53_only" degrade to "cluster_only" (knit-safe).
# Single source of truth for MMR-deficiency (MMRd) status, used by add_tcga_class() (the
# TCGA top tier), mmrd_aneuploidy_split() (report 08), and the report-09 heatmap MMR
# bar — so "MMRd" means the SAME thing in every report. Criterion: MMR-protein loss by IHC
# (MMR_class == deficient) ONLY. MSI-high (NGS) is deliberately NOT a criterion — it is a
# correlate. `mmr_class` is a character vector of MMR_class values; returns logical, NA
# where the MMR IHC status is unknown (so callers can decide how to treat unknowns).
is_mmrd <- function(mmr_class, mmr_deficient = attend_levels$mmr_deficient) {
  as.character(mmr_class) == mmr_deficient
}

add_tcga_class <- function(df,
                           msi_class_col = "MSI_class",   # kept for back-compat; NOT a criterion
                           mmr_class_col = "MMR_class",
                           cnv_high_col  = "cnv_high",
                           tp53_col      = "TP53_pathogenic",
                           serous_rule   = attend_tcga_serous_rule,
                           lev           = attend_tcga_levels,
                           mmr_lev       = attend_levels) {
  n  <- nrow(df)
  lg <- function(col) if (col %in% names(df)) as.logical(df[[col]]) else rep(NA, n)
  cnh  <- lg(cnv_high_col)

  # MMRd = MMR-protein loss by IHC ONLY (MMRd vs MMRp). MSI-high (NGS) is NOT used as a
  # criterion (a per-project decision so the top tier is driven purely by MMR IHC). A
  # sample with unknown MMR IHC (NA) is treated as NOT MMRd, so it falls through to the
  # CN-high/CN-low split rather than the top tier.
  mmr  <- if (mmr_class_col %in% names(df)) is_mmrd(df[[mmr_class_col]], mmr_lev$mmr_deficient) else rep(NA, n)
  mmrd <- !is.na(mmr) & mmr

  tp53      <- lg(tp53_col)
  have_tp53 <- tp53_col %in% names(df)
  if (!have_tp53 && serous_rule %in% c("and", "tp53_only")) {
    message("add_tcga_class(): '", tp53_col, "' absent — serous rule '", serous_rule,
            "' falls back to the CNV cluster alone.")
    serous_rule <- "cluster_only"
  }
  serous <- switch(serous_rule,
    and          = !is.na(cnh) & cnh & !is.na(tp53) & tp53,   # cluster AND TP53 (opt-in, higher precision)
    or           = (!is.na(cnh) & cnh) | (!is.na(tp53) & tp53),
    cluster_only = !is.na(cnh) & cnh,
    tp53_only    = !is.na(tp53) & tp53,
    stop("add_tcga_class(): unknown serous_rule '", serous_rule, "'"))

  cls <- rep(NA_character_, n)
  cls[!is.na(cnh)] <- lev$cn_low            # CNV-profiled but not serous -> endometrioid
  cls[serous]      <- lev$cn_high           # serous-like, per `serous_rule`
  cls[mmrd]        <- lev$msi               # MMRd is the top of the cascade (no POLE)
  df$tcga_class <- factor(cls, levels = c(lev$msi, lev$cn_low, lev$cn_high))
  df
}

# =============================================================================
# --- 13. Ancestry-aware TMB recompute (report 01) ---------------------------
# Tumor-only TMB is inflated by residual germline variants, and the inflation is
# ANCESTRY-DEPENDENT because population reference panels (gnomAD) are Euro-skewed
# (Nassar et al., Cancer Cell 2022). The ATTEND Funcotator MAFs carry gnomAD
# exome v4.1 PER-SUB-POPULATION allele frequencies, so we can recompute TMB in R
# with an ancestry-appropriate germline filter — no upstream re-run.
#
# The BASELINE (status-quo, no ancestry correction) is the upstream tmb__TMB_SCORE.
# Two ancestry-aware TMBs are recomputed on top of it, differing ONLY in the germline
# filter (both calibrated to the TMB_SCORE scale via the internal unfiltered count):
#   robust   — drop variants common in ANY population (grpmax >= threshold);
#              ancestry-robust, needs NO ancestry label
#   matched  — drop variants common in the PATIENT's own ancestry sub-population;
#              most faithful, needs a per-patient ancestry label
# (An unfiltered "raw" count is computed internally only to set the exome denominator
# — it is NOT emitted, since it duplicates tmb__TMB_SCORE.) A variant is germline
# (dropped) iff its chosen gnomAD AF >= threshold; rare/novel are kept as candidate somatic.
attend_tmb_recompute <- list(
  af_threshold    = 0.001,                 # gnomAD AF >= this => germline
  class_col       = "Variant_Classification",
  gene_col        = "Hugo_Symbol",
  grpmax_col      = "gnomAD_exome_AF_grpmax",
  subpop_cols     = c(afr = "gnomAD_exome_AF_afr", amr = "gnomAD_exome_AF_amr",
                      asj = "gnomAD_exome_AF_asj", eas = "gnomAD_exome_AF_eas",
                      fin = "gnomAD_exome_AF_fin", nfe = "gnomAD_exome_AF_nfe",
                      sas = "gnomAD_exome_AF_sas", mid = "gnomAD_exome_AF_mid"),
  counted_classes = c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
                      "Splice_Site", "Translation_Start_Site", "Frame_Shift_Del",
                      "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins"),
  default_exome_mb = 34                     # WES coding target (~Mb) if not calibrated
)

# The master columns the recomputed TMBs land in (tmb__ prefix convention).
# The BASELINE everywhere is the upstream `tmb__TMB_SCORE` (attend_cols$tmb) — the
# `raw` (unfiltered) recompute is NOT emitted as a column because it duplicates the
# upstream score by construction (it is calibrated to the same scale); it survives
# only as an internal count in tmb_variant_table() to set the exome denominator.
# `nassar` is the ancestry-specific regression recalibration (see §14).
# EXACTLY TWO TMB DEFINITIONS ARE EVER PLOTTED, IN EVERY REPORT:
#
#   normal             tmb__TMB_SCORE   the upstream tumour-only score, as delivered
#   ancestry corrected tmb__TMB_nassar  that score after the ancestry recalibration
#
# `tmb__TMB_robust` and `tmb__TMB_matched` are still COMPUTED in report 01 and still
# live on the master, because the ancestry correction is fitted on top of the robust
# recompute (attend_tmb_nassar$source_col). They are INTERMEDIATES: never plotted,
# never faceted over, never labelled. tmb_defs_present() cannot return them, so a
# report cannot reintroduce them by accident.
attend_cols$tmb_intermediate <- c(robust  = "tmb__TMB_robust",
                                  matched = "tmb__TMB_matched")
attend_cols$tmb_set <- c(nassar = "tmb__TMB_nassar")

# The reader-facing name of each definition. Reports label every facet, axis, and
# table row through this map rather than printing the internal name, so "nassar"
# never reaches the knitted page.
attend_tmb_labels <- c(upstream = "TMB (normal)",
                       nassar   = "TMB (ancestry corrected)")

# Look up display labels for the definitions tmb_defs_present() returned. Falls back
# to the column name for anything unmapped, so this can never error a knit.
tmb_def_labels <- function(defs, labels = attend_tmb_labels) {
  out <- unname(labels[names(defs)])
  ifelse(is.na(out), unname(defs), out)
}

# Syntactically-safe variants of the same names, for places that need the label as a
# COLUMN name rather than free text — the maftools oncoplot draws its clinicalFeatures
# column names directly onto the figure, so "TMB_class.nassar" would put the internal
# name on the page. Same map, no spaces or parentheses.
attend_tmb_track_names <- c(upstream = "TMB_normal",
                            nassar   = "TMB_ancestry_corrected")
tmb_def_track_names <- function(defs, names_map = attend_tmb_track_names) {
  out <- unname(names_map[names(defs)])
  ifelse(is.na(out), make.names(unname(defs)), out)
}

# --- TMB definition helpers (used by EVERY report that plots TMB) ------------
# One source of truth for "reproduce this figure under each TMB definition".
# The design intent is that no TMB plot hard-codes a single column: reports loop
# over `tmb_defs_present(df)` instead. This centralizes the idiom that reports
# 04/09 previously inlined so 03/06/07 cannot drift from it.
#
#   tmb_defs_present(df) -> ordered, NAMED vector of the TMB columns actually on
#     `df`: `tmb__TMB_SCORE` (normal) first, then `tmb__TMB_nassar` (ancestry
#     corrected) if present. Nassar is absent/all-NA until the ancestry
#     `clinical_col` is set, so callers that filter NA rows degrade gracefully.
#   tmb_class_of(x) -> binarise any TMB vector into "TMB-high"/"TMB-low" at the
#     configured cutoff, so TMB_class-keyed figures (KM, oncoplot) — which the
#     upstream-only `TMB_class` column cannot express — can be redrawn per
#     definition. Mirrors the derivation in add_molecular_classes().
tmb_defs_present <- function(df, cols = attend_cols) {
  defs <- c(upstream = cols$tmb, cols$tmb_set)
  defs[defs %in% names(df)]
}
tmb_class_of <- function(x, thr = attend_thresholds$tmb) {
  ifelse(suppressWarnings(as.numeric(x)) >= thr, "TMB-high", "TMB-low")
}

# Self-reported clinical ancestry -> gnomAD sub-population code (for MATCHED + NASSAR TMB).
# Matching is EXACT on the FULL lower-cased, trimmed label — ancestry_by_pid() does a
# whole-string lookup, NOT substring/token — so every map entry must be a COMPLETE label
# exactly as it appears in the data. Any label not listed here resolves to NA (that
# patient then gets no matched/nassar TMB, but everything else still runs; knit-safe).
#
# Set for the ATTEND cohort: `gianlu__RACE` has three values — White (210) / Asian (24)
# / Black Or African American (1). "Asian" is UNQUALIFIED in the source, so it is mapped
# to gnomAD EAS (East Asian) by default; change "asian" to the `sas` bucket below if this
# cohort's Asian patients are South Asian. Other labels are kept for future data.
attend_ancestry <- list(
  clinical_col = "gianlu__RACE",
  map = list(
    nfe = c("white", "caucasian", "european", "europe", "italian", "non-finnish european"),
    fin = c("finnish"),
    afr = c("black", "african", "african american", "africa", "black or african american"),
    eas = c("east asian", "chinese", "japanese", "korean", "asian"),
    sas = c("south asian", "indian", "pakistani", "bangladeshi"),
    amr = c("hispanic", "latino", "latina", "admixed american", "south american"),
    asj = c("ashkenazi", "jewish", "ashkenazi jewish")
  )
)

# Map a data frame with pid + a self-reported ancestry column to pid + subpop code.
# Returns NULL if the column is absent (matched TMB then stays NA — knit-safe).
ancestry_by_pid <- function(df, cfg = attend_ancestry) {
  col <- cfg$clinical_col
  if (is.null(col) || !col %in% names(df)) {
    message("ancestry_by_pid(): '", col, "' not in data — matched TMB unavailable."); return(NULL)
  }
  code_of <- unlist(lapply(names(cfg$map),
    function(code) setNames(rep(code, length(cfg$map[[code]])), tolower(cfg$map[[code]]))))
  df |>
    transmute(pid, .lab = tolower(trimws(as.character(.data[[col]])))) |>
    mutate(subpop = unname(code_of[.lab])) |>
    filter(!is.na(subpop)) |>
    distinct(pid, subpop)
}

# Recompute the three TMBs per patient from a gnomAD-annotated long MAF.
#   maf_tmb    : long MAF with Tumor_Sample_Barcode, class/gene cols, gnomAD AF cols
#                (load_maf_tmb() in load_wes_results.R)
#   crosswalk  : barcode -> pid (build_barcode_pid())
#   anc_by_pid : pid + subpop (ancestry_by_pid()); NULL -> matched = NA
#   anchor_tmb : optional pid + tmb_anchor (upstream tmb__TMB_SCORE) to CALIBRATE the
#                exome-size denominator so the recomputed counts sit on the upstream
#                scale (robust/matched then share a denominator; the only difference
#                from the baseline TMB_SCORE is the germline filter). Falls back to
#                default_exome_mb.
# Returns pid + tmb__TMB_robust / _matched (mut/Mb). NULL if no MAF.
tmb_variant_table <- function(maf_tmb, crosswalk, anc_by_pid = NULL,
                              anchor_tmb = NULL, cfg = attend_tmb_recompute) {
  if (is.null(maf_tmb) || nrow(maf_tmb) == 0 || is.null(crosswalk)) return(NULL)

  coding <- maf_tmb |>
    filter(.data[[cfg$class_col]] %in% cfg$counted_classes) |>
    mutate(ID = as.character(Tumor_Sample_Barcode))
  if (nrow(coding) == 0) return(NULL)

  keep_rare <- function(af) is.na(af) | af < cfg$af_threshold

  n_raw <- coding |> count(ID, name = "n_raw")

  n_rob <- if (cfg$grpmax_col %in% names(coding)) {
    coding |> mutate(.af = suppressWarnings(as.numeric(.data[[cfg$grpmax_col]]))) |>
      filter(keep_rare(.af)) |> count(ID, name = "n_robust")
  } else { n_raw |> transmute(ID, n_robust = n_raw) }

  n_mat <- NULL
  if (!is.null(anc_by_pid)) {
    anc_id <- anc_by_pid |> inner_join(crosswalk, by = "pid") |> transmute(ID = barcode, subpop)
    cm <- coding |> inner_join(anc_id, by = "ID")
    if (nrow(cm) > 0) {
      n_mat <- purrr::map_dfr(names(cfg$subpop_cols), function(code) {
        col <- cfg$subpop_cols[[code]]
        if (!col %in% names(cm)) return(tibble())
        cm |> filter(subpop == code) |>
          mutate(.af = suppressWarnings(as.numeric(.data[[col]]))) |>
          filter(keep_rare(.af))
      }) |> count(ID, name = "n_matched")
    }
  }

  counts <- n_raw |> full_join(n_rob, by = "ID")
  counts <- if (!is.null(n_mat)) full_join(counts, n_mat, by = "ID") else
              mutate(counts, n_matched = NA_real_)
  counts <- counts |> mutate(across(c(n_raw, n_robust), ~ replace_na(as.numeric(.x), 0)))

  perpid <- counts |>
    inner_join(crosswalk, by = c("ID" = "barcode")) |>
    group_by(pid) |>
    summarise(across(c(n_raw, n_robust, n_matched), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
    mutate(across(c(n_raw, n_robust, n_matched), ~ ifelse(is.nan(.x), NA_real_, .x)))

  # calibrate the exome-size denominator to the upstream TMB scale if given
  mb <- cfg$default_exome_mb
  if (!is.null(anchor_tmb) && "tmb_anchor" %in% names(anchor_tmb)) {
    a <- perpid |> inner_join(anchor_tmb, by = "pid")
    mr <- stats::median(a$n_raw, na.rm = TRUE); ma <- stats::median(a$tmb_anchor, na.rm = TRUE)
    if (is.finite(mr) && is.finite(ma) && ma > 0) mb <- mr / ma
  }

  # `raw` (= n_raw / mb) is NOT returned: it duplicates the upstream tmb__TMB_SCORE
  # (which is the baseline), and n_raw was only needed above to calibrate `mb`.
  perpid |> transmute(pid,
                      `tmb__TMB_robust`  = n_robust  / mb,
                      `tmb__TMB_matched` = n_matched / mb)
}

# =============================================================================
# --- 14. Nassar ancestry regression recalibration (report 00-methods) ---------------
# The headline method of Nassar et al. (Cancer Cell 2022) is NOT a germline
# filter — it is an ancestry-specific AFFINE recalibration of the TMB value:
#   TMB_corrected = m_group * TMB_tumor-only + b_group,
# with (m, b) fit separately for European vs non-European against matched
# tumor-normal (paired) TMB. Here it is applied on top of the grpmax-filtered
# `robust` TMB (which equals their maxPOPAF <= 0.1% germline filter).
#
# >>> IMPORTANT — recalibration DEFAULTS TO IDENTITY (active_source = "identity") <<<
# so `tmb__TMB_nassar` == `tmb__TMB_robust` until you pick a coefficient set.
# The paper's PUBLISHED numbers (in `paper_coef`) are ONCOPANEL-specific and do NOT
# strictly transfer to whole-exome data — they are provided as selectable presets
# (`active_source = "oncopanel_v3"` etc.) for anyone who explicitly wants them, but
# are OFF by default. The faithful path for ATTEND (WES) is to FIT (m, b) per ancestry
# on TCGA (paired-WES TMB ~ tumor-only-WES TMB with the same robust filter) via
# fit_nassar_coefficients(), paste into `coef`, and set active_source = "custom".
attend_tmb_nassar <- list(
  # gnomAD sub-populations counted as "European" for the binary Nassar split
  # (nfe = non-Finnish European, fin = Finnish, asj = Ashkenazi; edit as needed).
  european_subpops = c("nfe", "fin", "asj"),
  source_col       = "tmb__TMB_robust",     # recalibrate on top of the robust TMB

  # Which coefficient set apply_nassar_recalibration() applies (via nassar_active_coef()):
  #   "identity"      -> m=1, b=0  => tmb__TMB_nassar == robust (SAFE DEFAULT)
  #   "oncopanel_v3" / "_v2" / "_v1" -> the PUBLISHED Nassar et al. 2022 coefficients in
  #        `paper_coef`. ⚠️ Fit on OncoPanel (targeted-panel) tumor-only TMB, NOT WES —
  #        applying them to ATTEND's exome robust TMB is an off-scale approximation (the
  #        slope maps OncoPanel->WES-truth, not WES->WES-truth). Provided because it was
  #        requested; use knowingly. Set to "oncopanel_v3" to activate the paper's values.
  #   "custom"        -> whatever is in `coef` below (e.g. your own WES-fit values)
  # Set to the paper's most-recent panel per request; revert to "identity" for nassar==robust,
  # or "custom" once you have WES-fit coefficients.
  active_source = "oncopanel_v3",

  coef = list(                              # used when active_source == "custom"
    european     = list(m = 1, b = 0),
    non_european = list(m = 1, b = 0)
  ),

  # PUBLISHED coefficients — Nassar et al., Cancer Cell 2022 (STAR Methods, "TMB
  # Recalibration"): TMB_paired ~ m * TMB_tumor_only + b, fit per ancestry per OncoPanel
  # version with TCGA WES paired tumor/normal as ground truth. European vs Non-European.
  paper_coef = list(
    oncopanel_v1 = list(european = list(m = 0.989, b = -2.18), non_european = list(m = 0.821, b = -1.71)),
    oncopanel_v2 = list(european = list(m = 0.994, b = -2.12), non_european = list(m = 0.840, b = -1.71)),
    oncopanel_v3 = list(european = list(m = 1.094, b = -1.94), non_european = list(m = 0.895, b = -1.29))
  ),

  # Optional reference cohort to FIT coefficients from (report 00-methods "own coefficients"
  # section). Drop a table with paired (tumor-normal) + tumor-only TMB + an ancestry
  # `group` column ("european"/"non_european") into data/nassar/. Absent -> the report
  # just reports the active coefficients and explains how to fit. Knit-safe.
  reference = list(
    dir        = "nassar",                  # data/nassar/
    glob       = "*reference*",             # any *reference*.tsv/.csv
    paired     = "paired_tmb",
    tumor_only = "tumor_only_tmb",
    group      = "group"
  )
)

# Resolve the coefficient set apply_nassar_recalibration() should use, per active_source
# ("identity" | "oncopanel_v1/2/3" | "custom"). Unknown source -> identity (knit-safe).
nassar_active_coef <- function(cfg = attend_tmb_nassar) {
  src <- if (is.null(cfg$active_source)) "identity" else cfg$active_source
  identity_coef <- list(european = list(m = 1, b = 0), non_european = list(m = 1, b = 0))
  if (src == "identity") return(identity_coef)
  if (src == "custom")   return(cfg$coef)
  if (!is.null(cfg$paper_coef[[src]])) return(cfg$paper_coef[[src]])
  message("nassar_active_coef(): unknown active_source '", src, "' — using identity.")
  identity_coef
}

# Tidy the ACTIVE Nassar coefficients into a table + an `identity` flag and the source
# name (report 00-methods). `identity` TRUE => tmb__TMB_nassar == robust; that gates the
# recalibration-effect panels (they only render once non-identity coefficients apply).
nassar_coef_table <- function(cfg = attend_tmb_nassar) {
  co <- nassar_active_coef(cfg)
  gs <- names(co)
  tab <- tibble::tibble(
    group = gs,
    m     = vapply(gs, function(g) co[[g]]$m, numeric(1)),
    b     = vapply(gs, function(g) co[[g]]$b, numeric(1))
  )
  list(table    = tab,
       identity = all(tab$m == 1 & tab$b == 0),
       source   = if (is.null(cfg$active_source)) "identity" else cfg$active_source)
}

# Add `tmb__TMB_nassar` = m_group * source + b_group, grouped European vs
# non-European from a pid->subpop table (ancestry_by_pid()). Patients without an
# ancestry label get NA. Knit-safe: missing source column -> unchanged data frame.
apply_nassar_recalibration <- function(df, anc_by_pid = NULL, cfg = attend_tmb_nassar) {
  src <- cfg$source_col
  if (!src %in% names(df)) {
    message("apply_nassar_recalibration(): '", src, "' absent — skipping."); return(df)
  }
  if (is.null(anc_by_pid) || nrow(anc_by_pid) == 0) {
    df$`tmb__TMB_nassar` <- NA_real_; return(df)
  }
  grp <- anc_by_pid |>
    mutate(group = ifelse(subpop %in% cfg$european_subpops, "european", "non_european")) |>
    distinct(pid, group) |>
    group_by(pid) |> slice(1) |> ungroup()          # one group per patient

  co <- nassar_active_coef(cfg)             # identity / oncopanel_vX / custom
  d <- df |> left_join(grp, by = "pid")
  m <- ifelse(d$group == "european",     co$european$m,
       ifelse(d$group == "non_european", co$non_european$m, NA_real_))
  b <- ifelse(d$group == "european",     co$european$b,
       ifelse(d$group == "non_european", co$non_european$b, NA_real_))
  df$`tmb__TMB_nassar` <- m * suppressWarnings(as.numeric(d[[src]])) + b
  df
}

# Fit the Nassar coefficients from a reference cohort that has BOTH paired
# (tumor-normal) and tumor-only TMB plus an ancestry group — i.e. reproduce the
# paper's regression for YOUR platform. Run on TCGA WES: compute paired-WES TMB and
# tumor-only-WES TMB (tumor-only simulated with the same robust/grpmax 0.1% filter),
# label European vs non-European, then:
#     coefs <- fit_nassar_coefficients(tcga_ref)
# and paste `coefs` into attend_tmb_nassar$coef. Returns a list with $european /
# $non_european, each list(m, b); groups with < 3 samples fall back to identity.
fit_nassar_coefficients <- function(ref, paired = "paired_tmb",
                                    tumor_only = "tumor_only_tmb", group = "group") {
  fit_one <- function(d) {
    d <- d[stats::complete.cases(d[c(paired, tumor_only)]), , drop = FALSE]
    if (nrow(d) < 3) return(list(m = 1, b = 0))
    co <- stats::coef(stats::lm(d[[paired]] ~ d[[tumor_only]]))
    list(m = unname(co[2]), b = unname(co[1]))
  }
  out <- lapply(split(ref, ref[[group]]), fit_one)
  # normalise names to what apply_nassar_recalibration() expects
  if (!"european" %in% names(out))     out$european     <- list(m = 1, b = 0)
  if (!"non_european" %in% names(out)) out$non_european <- list(m = 1, b = 0)
  out[c("european", "non_european")]
}

# =============================================================================
# --- 15. MMRd aneuploidy subclassification --------------------------------
# NOTE: currently UNUSED by any report. mmrd_aneuploidy_split() lost its caller when
# the MMRd-subclass and cross-cohort reports were removed; dip_bimodality() below is
# still used by report 07. Kept because both are tested and cheap to re-wire.
# The CUSTOM MMRd split (MMRd_subclassification_workflow §1.2): the custom feature is
# the aneuploidy score, and we divide the MMR-deficient group into aneuploidy-HIGH vs
# aneuploidy-LOW. MMRd membership = MSI-high (NGS) OR MMR-protein loss (IHC) — the same
# definition the cascade uses (add_tcga_class()). The cut is chosen principledly: the
# pre-registered attend_thresholds$aneuploidy by default (a frozen quantile-style rule,
# §1.2's endorsed portable option), optionally refined to the 2-component Gaussian-
# mixture boundary when `mclust` is present and the score is bimodal. Hartigan's dip
# test (`diptest`) quantifies whether a split even EXISTS (bimodal) vs a cut through
# unimodal noise. Returns the MMRd rows with:
#   aneu_value     numeric aneuploidy score used
#   aneu_subclass  factor "MMRd aneuploidy-low" / "MMRd aneuploidy-high"
# and attributes: cut, method ("fixed"/"gmm"), dip_p, n_mmrd. Knit-safe (all-NA in ->
# 0-row out). NOTE (workflow §5, circularity): because the split is DEFINED by
# aneuploidy, do NOT treat aneuploidy, SCNA burden, HRD, or WGD as evidence the split
# is meaningful — those separate by construction. Meaningfulness must come from
# features NOT split on (survival, immune/IHC, specific mutations, an independent modality).
mmrd_aneuploidy_split <- function(df,
                                  aneu_col  = attend_cols$aneuploidy,
                                  msi_col   = "MSI_class",
                                  mmr_col   = "MMR_class",
                                  fixed_cut = attend_thresholds$aneuploidy,
                                  refine_gmm = TRUE,
                                  lev       = attend_levels) {
  n   <- nrow(df)
  # MMRd = MMR-protein loss by IHC ONLY (is_mmrd); MSI-high (NGS) is a correlate, not a
  # criterion — the single MMRd definition shared with add_tcga_class() so reports agree.
  mmr <- if (mmr_col %in% names(df)) is_mmrd(df[[mmr_col]], lev$mmr_deficient)   else rep(NA, n)
  mmrd <- !is.na(mmr) & mmr

  d <- df[mmrd, , drop = FALSE]
  a <- suppressWarnings(as.numeric(if (aneu_col %in% names(d)) d[[aneu_col]] else rep(NA, nrow(d))))
  d$aneu_value <- a
  ok <- is.finite(a)

  dip_p <- if (requireNamespace("diptest", quietly = TRUE) && sum(ok) >= 4)
    tryCatch(diptest::dip.test(a[ok])$p.value, error = function(e) NA_real_) else NA_real_

  cut <- fixed_cut; method <- "fixed"
  if (isTRUE(refine_gmm) && requireNamespace("mclust", quietly = TRUE) && sum(ok) >= 6) {
    mc <- tryCatch(mclust::Mclust(a[ok], G = 2, verbose = FALSE), error = function(e) NULL)
    if (!is.null(mc) && length(unique(mc$classification)) == 2) {
      lo_hi   <- order(mc$parameters$mean)                 # component means, low then high
      hi_comp <- lo_hi[2]
      is_hi   <- mc$classification == hi_comp
      # display boundary = midpoint of the gap between the two components
      cut     <- mean(c(max(a[ok][!is_hi], na.rm = TRUE), min(a[ok][is_hi], na.rm = TRUE)))
      if (is.finite(cut)) method <- "gmm" else cut <- fixed_cut
    }
  }

  d$aneu_subclass <- factor(ifelse(a >= cut, "MMRd aneuploidy-high", "MMRd aneuploidy-low"),
                            levels = c("MMRd aneuploidy-low", "MMRd aneuploidy-high"))
  attr(d, "cut")    <- cut
  attr(d, "method") <- method
  attr(d, "dip_p")  <- dip_p
  attr(d, "n_mmrd") <- sum(mmrd)
  d
}

# Bimodality of a numeric vector (Hartigan's dip test) — the scale-invariant "is there a
# split?" statistic for a cross-cohort contrast: a metastatic-cohort
# MMRd aneuploidy vector may be bimodal (a real split) while a primary-cohort one is not.
# Returns list(dip, p, n, bimodal). NA when `diptest` is absent or n < 4. PURE, no I/O.
dip_bimodality <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
  if (length(x) < 4 || !requireNamespace("diptest", quietly = TRUE))
    return(list(dip = NA_real_, p = NA_real_, n = length(x), bimodal = NA))
  dt <- tryCatch(diptest::dip.test(x), error = function(e) NULL)
  if (is.null(dt)) return(list(dip = NA_real_, p = NA_real_, n = length(x), bimodal = NA))
  list(dip = unname(dt$statistic), p = dt$p.value, n = length(x), bimodal = dt$p.value < 0.05)
}

# --- 16. TCGA UCEC reference cohort (report 07, TMB comparison) ---------------
# Cohort B for the MMRd aneuploidy contrast: TCGA UCEC (PRIMARY tumours) vs ATTEND
# (METASTATIC). Reads a cBioPortal PanCancer-Atlas clinical export dropped into
# data/tcga/ (study `ucec_tcga_pan_can_atlas_2018`; download via code/fetch_tcga_ucec.R
# or the study's bulk files data_clinical_{sample,patient}.txt). The cBioPortal clinical
# files carry 4 leading #-comment lines before the real header — the loader skips them.
# Column names below are the cBioPortal defaults; edit if your export differs.
attend_tcga_ref <- list(
  dir            = "tcga",                         # data/tcga/
  sample_glob    = "*clinical_sample*.txt",
  patient_glob   = "*clinical_patient*.txt",
  id_col         = "PATIENT_ID",
  aneuploidy_col = "ANEUPLOIDY_SCORE",             # Taylor et al. 2018 (0-39)
  subtype_col    = "SUBTYPE",                       # UCEC_MSI / _POLE / _CN_HIGH / _CN_LOW
  msi_score_col  = "MSI_SCORE_MANTIS",              # MANTIS; MSI-high ~ > 0.4
  msi_hi_mantis  = 0.4,
  fga_col        = "FRACTION_GENOME_ALTERED",
  tmb_col        = "TMB_NONSYNONYMOUS",
  os_time_col    = "OS_MONTHS",
  os_event_col   = "OS_STATUS",                     # "1:DECEASED"/"0:LIVING"
  msi_subtype_pattern = "MSI"                       # SUBTYPE marking MMRd/MSI tumours
)

# --- TCGA UCEC 2013 PUBLICATION cohort (report 07, TMB comparison) ------------
# The ORIGINAL Nature 2013 study (cBioPortal `ucec_tcga_pub`), NOT the 2018 PanCancer-Atlas
# restatement above. Two things make it the right source for reproducing Fig. 1a:
#   SUBTYPE        the paper's four integrated groups, spelled out in full
#                  ("POLE (Ultra-mutated)" / "MSI (Hyper-mutated)" /
#                   "Copy-number low (Endometriod)" [sic, TCGA's spelling] /
#                   "Copy-number high (Serous-like)") — patient-level, n = 232
#   CNA_CLUSTER_K4 the PUBLISHED Fig-1a copy-number clusters 1-4 — sample-level, n = 363.
#                  Cluster 4 is 60/60 concordant with "Copy-number high (Serous-like)",
#                  which is the empirical confirmation of attend_tcga_serous_rule =
#                  "cluster_only": TCGA's CN-high call IS cluster-4 membership, with no
#                  TP53 gate anywhere in it.
# The study has NO aneuploidy score (it predates Taylor et al. 2018); the fetcher joins
# one in from the PanCancer-Atlas study, hence `aneuploidy_file`. Genome build hg19.
attend_tcga_ref_2013 <- list(
  dir             = "tcga_2013",                    # data/tcga_2013/
  build           = "hg19",                         # SNP6 segments are hg19
  sample_glob     = "*clinical_sample*.txt",
  patient_glob    = "*clinical_patient*.txt",
  aneuploidy_file = "pancan_aneuploidy.tsv",        # PATIENT_ID + ANEUPLOIDY_SCORE
  seg_glob        = "*.seg",
  id_col          = "PATIENT_ID",
  sample_id_col   = "SAMPLE_ID",
  subtype_col     = "SUBTYPE",                      # the four 2013 integrated groups
  cluster_col     = "CNA_CLUSTER_K4",               # published Fig-1a clusters (1-4)
  # TCGA's OWN core/non-core flag, and the explanation for every missing SUBTYPE. Measured on
  # this export: DATA_CORE_SAMPLE is Y for 232 and N for 141 of 373, and the 141 non-core are
  # EXACTLY the 141 with no integrated subtype (0 disagreements either way). The integrated
  # classification needed the full multiplatform data, which only the core set has; the other
  # 141 carry SNP6 copy number and were never integratively classified. So "not classified" is
  # TCGA's study design, not missingness in our copy, and this column is how a report restricts
  # to the classified cohort without inventing a criterion.
  core_col        = "DATA_CORE_SAMPLE",             # Y/N — the 232-sample multiplatform core
  msi_call_col    = "MSI_STATUS_7_MARKER_CALL",     # MSS / MSI-H / MSI-L / Indeterminant
  mlh1_col        = "MLH1_SILENCING",               # 0/1
  mutrate_col     = "MUTATION_RATE_CLUSTER",        # 1_LOW / 2_HIGH / 3_HIGHEST
  grade_col       = "GRADE",                        # "Grade 1/2/3"  — complete (373/373)
  stage_col       = "TUMOR_STAGE_2009",             # "Stage I..IV" + "Unknown" (373/373)
  grade_levels    = c("Grade 1", "Grade 2", "Grade 3"),
  stage_levels    = c("Stage I", "Stage II", "Stage III", "Stage IV"),
  stage_unknown   = "Unknown",                      # mapped to NA, not to a 5th level
  fga_col         = "FRACTION_GENOME_ALTERED",      # study-native SCNA burden (n = 365)
  tmb_col         = "TMB_NONSYNONYMOUS",
  aneuploidy_col  = "ANEUPLOIDY_SCORE",             # retrofitted from PanCancer-Atlas
  # ASCETS run on the TCGA segments by code/run_ascets_tcga.sh — the PREFERRED aneuploidy
  # source, because it is the SAME statistic ATTEND uses (fraction of
  # evaluable arms altered, 0-1) rather than the Taylor et al. 2018 arm COUNT (0-39) that
  # cBioPortal serves. Same output file names as the ATTEND ASCETS run (attend_cnv$ascets),
  # but no id_strip — TCGA sample barcodes need no normalisation.
  ascets = list(
    dir             = "tcga_2013/ascets",           # data/tcga_2013/ascets/
    aneuploidy_glob = "*aneuploidy_scores*.txt",
    calls_glob      = "*arm_level_calls*.txt",
    armmeans_glob   = "*arm_weighted_average_segmeans*.txt"
  ),
  os_time_col     = "OS_MONTHS",
  os_event_col    = "OS_STATUS",
  dfs_time_col    = "DFS_MONTHS",                   # Fig-1b is progression-free survival
  dfs_event_col   = "DFS_STATUS",
  msi_subtype_pattern  = "MSI",
  # HOW MMRd IS CALLED. Default is TCGA-EXACT.
  #   "tcga_exact" (DEFAULT) — MMRd == the paper's own MSI (Hyper-mutated) integrative
  #        cluster, nothing else: n = 65. Samples with no SUBTYPE are NA (unclassified),
  #        NOT proficient — TCGA never called them either way. Use this whenever the claim
  #        is "TCGA's MSI group", i.e. anything quoted against the publication.
  #   "broad" — additionally counts a 7-marker MSI-H call in patients TCGA never
  #        classified (no exome -> no integrative cluster): n = 124. Trades fidelity for
  #        coverage; only defensible when the comparison needs MMR status across the whole
  #        cohort and the deviation is stated.
  mmrd_rule            = "tcga_exact",
  # POLE OUTRANKS MSI in the paper's Fig-2b cascade (POLE -> MSI -> CN-low/CN-high). Three
  # POLE-ultramutated tumours in this cohort are also MSI-H by the 7-marker assay, because
  # ultramutation destabilises microsatellites as a CONSEQUENCE of the mutation burden, not
  # through MMR loss. Without this pattern they would be called MMRd, contradicting the very
  # SUBTYPE column the same reports display.
  pole_subtype_pattern = "POLE",
  gistic_dir      = "tcga_2013/gistic/all",         # code/run_gistic_tcga.sh output
  # TCGA's OWN published GISTIC peaks — Supplementary Data File S2.1 of nature12113,
  # "all tumor amplifications" + "all tumor deletions" sheets, converted to TSV. 79 peaks
  # (43 amp / 36 del) over 5,046 genes. This is an EXTERNAL, PRE-REGISTERED feature set:
  # it was not fitted to ATTEND, so using it removes the cohort-fitted peak instability that
  # forces report 06 into leave-one-out stability testing, and it lets both cohorts be
  # clustered on identical features. The paper never published its GISTIC PARAMETERS
  # (Suppl. Methods says only "GISTIC 2.0"), so the peaks themselves are the only faithful
  # route to its feature space.
  published_peaks = "tcga_gistic_peaks_2013.tsv",
  # Level order for every plot/table: the cascade order of the paper's Fig. 2b.
  subtype_levels  = c("POLE (Ultra-mutated)", "MSI (Hyper-mutated)",
                      "Copy-number low (Endometriod)", "Copy-number high (Serous-like)"),
  # Short labels for axes where the full names do not fit.
  subtype_short   = c("POLE (Ultra-mutated)"           = "POLE",
                      "MSI (Hyper-mutated)"            = "MSI",
                      "Copy-number low (Endometriod)"  = "CN-low",
                      "Copy-number high (Serous-like)" = "CN-high")
)

# A `cnv`-SHAPED view of the TCGA GISTIC run, so the EXISTING GISTIC loaders work on it
# unchanged: load_gistic_thresholded(cnv = attend_tcga_gistic) and find_gistic_files() both
# resolve their folder from cnv$gistic$dir, and load_gistic_thresholded_at() strips ids with
# cnv$seg$id_strip. Nothing new to write — only the paths differ. The peak globs are copied
# from attend_cnv$gistic because GISTIC's output filenames are identical for any run.
attend_tcga_gistic <- list(
  gistic = list(
    dir                  = attend_tcga_ref_2013$gistic_dir,
    all_lesions_glob     = "*all_lesions*.txt",
    amp_genes_glob       = "*amp_genes*.txt",
    del_genes_glob       = "*del_genes*.txt",
    scores_glob          = "*scores.gistic",
    all_thresholded_glob = "*all_thresholded*.txt"
  ),
  # TCGA sample barcodes need no normalisation; "$^" is a regex that can never match, so
  # the loader's str_remove() is a no-op rather than silently mangling an id.
  seg = list(id_strip = "$^")
)

# =============================================================================
# §15  RECURRENT SCNA BY ANEUPLOIDY x MMR  (report 06)
#
# The 2x2 that report 06 tests. MMRd-high is the small cell (n=9 as of 2026-07-20),
# which is why the primary endpoint is an AGGREGATED panel score rather than a
# per-peak test — see the spec's power table. Panel loci are DIRECTIONAL: a deletion
# at MYC must not score as serous-like.
# =============================================================================
attend_scna <- list(

  # Pre-specified TCGA UCEC (Kandoth 2013) recurrent focal loci — Family A,
  # the confirmatory family. 8 amplifications, 4 deletions.
  panel = data.frame(
    locus     = c("MYC", "CCNE1", "ERBB2", "PIK3CA_SOX2", "SOX17",
                  "1q", "20q13", "FGFR3",
                  "PTEN", "RB1", "CDKN2A", "TP53"),
    cytoband  = c("8q24", "19q12", "17q12", "3q26", "8q11",
                  "1q", "20q13", "4p16",
                  "10q23", "13q14", "9p21", "17p13"),
    direction = c(rep("amp", 8), rep("del", 4)),
    stringsAsFactors = FALSE
  ),

  group_levels = c("MMRp-low", "MMRp-high", "MMRd-low", "MMRd-high"),

  maxseg             = 46000,   # hypersegmentation QC filter; exclusions are audited
  loo_min_overlap_bp = 1,       # leave-one-out peak "retained" = wide-limit overlap
  perm_B             = 10000,   # label-permutation replicates
  fdr_method         = "BH",    # Family B (discovery), genome-wide
  panel_p_adjust     = "holm",  # Family A (confirmatory), 12 loci

  # Data-driven panel size — Family B ONLY. The loci are chosen from the group
  # labels, so this panel can never enter Family A no matter how it scores; see
  # perm_test_selected_panel(), which nests the selection inside the permutation.
  # Matched to the 12 pre-specified loci so the two scores are on one scale.
  select_n           = 12,
  # Replicates for the NESTED test. Lower than perm_B because every replicate
  # re-runs the selection: cost is B x (selection + scoring), not B x scoring.
  select_perm_B      = 2000
)

#' Cross MMR status with aneuploidy class into the report-15 grouping factor.
#'
#' Composes existing derived columns rather than re-deriving either one: MMR status
#' comes from is_mmrd() (IHC only, MSI deliberately excluded) and aneuploidy_class
#' from add_molecular_classes(). Re-deriving would let report 06 drift from 09/11/12.
#'
#' NA MMR status propagates to NA in both output columns.
add_scna_group <- function(df,
                           mmr_col        = "MMR_class",
                           aneu_class_col = "aneuploidy_class",
                           lev            = attend_levels,
                           cfg            = attend_scna) {
  if (!all(c(mmr_col, aneu_class_col) %in% names(df))) {
    stop("add_scna_group(): missing column(s): ",
         paste(setdiff(c(mmr_col, aneu_class_col), names(df)), collapse = ", "))
  }

  d <- is_mmrd(df[[mmr_col]], lev$mmr_deficient)          # TRUE / FALSE / NA
  mmr <- ifelse(is.na(d), NA_character_, ifelse(d, "MMRd", "MMRp"))

  hl <- ifelse(grepl("high", as.character(df[[aneu_class_col]]), fixed = TRUE),
               "high", "low")
  hl[is.na(df[[aneu_class_col]])] <- NA_character_

  df$mmr_group  <- factor(mmr, levels = c("MMRp", "MMRd"))
  df$scna_group <- factor(ifelse(is.na(mmr) | is.na(hl), NA_character_,
                                 paste0(mmr, "-", hl)),
                          levels = cfg$group_levels)
  df
}
