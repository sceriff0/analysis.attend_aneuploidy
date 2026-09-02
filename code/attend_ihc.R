# =============================================================================
# attend_ihc.R  —  FlowPath IHC derived metrics (used by report 02, 30 and 32).
#
# The "spatial/IHC" readout is built from the per-image FlowPath cell counts
# (load_ihc_data()) joined to the imaging pathologist annotations
# (load_imaging_data(), which carries image_id + ID = pid). Defining the ratio
# metrics here means report 02 (clinical<->IHC concordance) and report 05
# (IHC vs aneuploidy) cannot drift apart.
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))

# Per-image: join FlowPath counts to imaging annotations (by image_id) and derive
# positivity / composition ratios. NA propagates where a batch lacks PD-1/PD-L1.
ihc_image_metrics <- function(ihc_data, imaging_data) {
  ihc_data |>
    dplyr::left_join(imaging_data, by = "image_id") |>
    dplyr::mutate(
      tumor_over_all                 = n_tumor_total        / n_cells,
      cd45_over_inside               = n_cd45_inside        / n_inside,
      tumor_aridia_over_tumor_inside = n_tumor_aridia_inside / n_tumor_inside,
      tumor_pdl1_over_tumor_inside   = n_tumor_pdl1_inside  / n_tumor_inside,
      cd45_pd1_over_cd45_inside      = n_cd45_pd1_inside    / n_cd45_inside,
      cd45_pdl1_over_cd45_inside     = n_cd45_pdl1_inside   / n_cd45_inside
    )
}

# The FlowPath ratio metric columns (one value per image, then per patient).
ihc_metric_cols <- c(
  "tumor_over_all", "cd45_over_inside", "tumor_aridia_over_tumor_inside",
  "tumor_pdl1_over_tumor_inside", "cd45_pd1_over_cd45_inside", "cd45_pdl1_over_cd45_inside"
)

# Collapse to one row per patient (mean across that patient's images). imaging's
# `ID` column is the Unique Subject Identifier == pid.
ihc_patient_metrics <- function(ihc_data, imaging_data, metric_cols = ihc_metric_cols) {
  ihc_image_metrics(ihc_data, imaging_data) |>
    dplyr::filter(!is.na(ID)) |>
    dplyr::group_by(pid = as.character(ID)) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(metric_cols),
                                   ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

# Per-patient cell-type composition inside the annotation, normalised three ways:
# by all cells inside, by tumour cells inside, and by CD45+ cells inside. Pools a
# patient's images (sum counts, sum denominators). `ihc_celltypes` comes from
# load_ihc_celltypes(); imaging maps image_id -> pid (ID). Long: pid, cell_type,
# n_cell, frac_inside, frac_tumor, frac_cd45.
ihc_celltype_metrics <- function(ihc_celltypes, imaging_data) {
  img <- imaging_data |>
    dplyr::distinct(image_id, pid = as.character(ID)) |>
    dplyr::filter(!is.na(pid))
  ct <- ihc_celltypes |> dplyr::inner_join(img, by = "image_id")

  cells <- ct |>
    dplyr::group_by(pid, cell_type) |>
    dplyr::summarise(n_cell = sum(n_cell), .groups = "drop")
  denom <- ct |>
    dplyr::distinct(pid, image_id, n_inside, n_tumor_inside, n_cd45_inside) |>
    dplyr::group_by(pid) |>
    dplyr::summarise(dplyr::across(c(n_inside, n_tumor_inside, n_cd45_inside), sum),
                     .groups = "drop")

  cells |>
    dplyr::left_join(denom, by = "pid") |>
    dplyr::mutate(frac_inside = n_cell / n_inside,
                  frac_tumor  = n_cell / n_tumor_inside,
                  frac_cd45   = n_cell / n_cd45_inside)
}

# TISSUE CONTENT: how much of the annotation is tumour, and how much is leukocyte, on one
# axis. The two do not come from one table -- `Tumor` is a cell_type row in
# ihc_celltype_metrics(), while CD45+ is a per-image DENOMINATOR constant (n_cd45_inside)
# that only ihc_immune_metrics() surfaces -- so this unions them onto a shared `cell_type`
# column and returns pid / cell_type / frac, the shape the composition box helpers already
# facet on.
#
# `all cells inside` is the only denominator where BOTH are meaningful, which is why it is
# the default rather than a loop over all three: Tumor/tumour_inside is identically 1, and
# CD45+/cd45_inside is identically 1 and is already dropped upstream. A caller asking for
# another denominator gets whichever series survives, and a message naming the one that did
# not -- silently drawing a panel of 1.0 would look like a result.
#
# Returns a zero-row frame (not NULL) when neither series is present, so a call site can
# nrow()-check it the same way it checks every other composition table.
ihc_content_metrics <- function(celltype_metrics, immune_metrics,
                                denom_label = "all cells inside",
                                frac_col    = "frac_inside",
                                tumor_cell_type = "Tumor",
                                immune_metric   = "CD45+") {
  empty <- data.frame(pid = character(0), cell_type = character(0), frac = numeric(0))

  tum <- if (!is.null(celltype_metrics) && nrow(celltype_metrics) &&
             frac_col %in% names(celltype_metrics)) {
    d <- celltype_metrics[celltype_metrics$cell_type == tumor_cell_type, , drop = FALSE]
    if (nrow(d)) data.frame(pid = d$pid, cell_type = tumor_cell_type,
                            frac = d[[frac_col]], stringsAsFactors = FALSE) else empty
  } else empty

  imm <- if (!is.null(immune_metrics) && nrow(immune_metrics)) {
    d <- immune_metrics[immune_metrics$metric == immune_metric &
                        immune_metrics$denom_label == denom_label, , drop = FALSE]
    if (nrow(d)) data.frame(pid = d$pid, cell_type = immune_metric,
                            frac = d$frac, stringsAsFactors = FALSE) else empty
  } else empty

  if (!nrow(tum)) message("Content panel: no '", tumor_cell_type, "' rows for ", denom_label, ".")
  if (!nrow(imm)) message("Content panel: no '", immune_metric, "' rows for ", denom_label, ".")

  out <- rbind(tum, imm)
  # Tumour first, then immune — the reader's order, and it must not depend on rbind().
  out$cell_type <- factor(out$cell_type, levels = c(tumor_cell_type, immune_metric))
  out[is.finite(out$frac), , drop = FALSE]
}

# Arcsine-square-root transform for a proportion in [0, 1] — the classic
# variance-stabiliser for fraction data, which is what makes a parametric t-test
# defensible on these composition fractions. Clamps out-of-range inputs.
arcsin_sqrt <- function(p) asin(sqrt(pmin(pmax(p, 0), 1)))

# Per-patient IMMUNE content inside the annotation, for report 05: the two
# leukocyte metrics (CD45+ total and CD3+CD45+ double-positive T cells) under each
# normalisation, compared later between aneuploidy classes. Pools a patient's
# images (sum counts, sum denominators), same as ihc_celltype_metrics(). Reads the
# per-image constants (n_inside/n_tumor_inside/n_cd45_inside/n_cd3cd45_inside) that
# load_ihc_celltypes() carries. The trivial CD45+/CD45+ self-normalisation (==1) is
# dropped. Long: pid, metric, denom_label, n_pos, frac. When CD3 staining is absent
# (n_cd3cd45_inside all NA) the CD3+CD45+ rows come through NA and downstream skips.
ihc_immune_metrics <- function(ihc_celltypes, imaging_data) {
  img <- imaging_data |>
    dplyr::distinct(image_id, pid = as.character(ID)) |>
    dplyr::filter(!is.na(pid))

  # per-image constants, then pooled per patient
  per_pt <- ihc_celltypes |>
    dplyr::inner_join(img, by = "image_id") |>
    dplyr::distinct(pid, image_id, n_inside, n_tumor_inside, n_cd45_inside,
                    n_cd3cd45_inside) |>
    dplyr::group_by(pid) |>
    dplyr::summarise(dplyr::across(c(n_inside, n_tumor_inside, n_cd45_inside,
                                     n_cd3cd45_inside), ~ sum(.x, na.rm = FALSE)),
                     .groups = "drop")

  # (metric, count, denom column, denom label) — CD45+/CD45+ intentionally omitted
  specs <- tibble::tribble(
    ~metric,       ~num,               ~den,             ~denom_label,
    "CD45+",       "n_cd45_inside",    "n_inside",       "all cells inside",
    "CD45+",       "n_cd45_inside",    "n_tumor_inside", "tumour cells inside",
    "CD3+CD45+",   "n_cd3cd45_inside", "n_inside",       "all cells inside",
    "CD3+CD45+",   "n_cd3cd45_inside", "n_tumor_inside", "tumour cells inside",
    "CD3+CD45+",   "n_cd3cd45_inside", "n_cd45_inside",  "CD45+ cells inside"
  )

  purrr::pmap_dfr(specs, function(metric, num, den, denom_label) {
    per_pt |>
      dplyr::transmute(pid, metric = metric, denom_label = denom_label,
                       n_pos = .data[[num]],
                       frac  = .data[[num]] / .data[[den]])
  })
}
