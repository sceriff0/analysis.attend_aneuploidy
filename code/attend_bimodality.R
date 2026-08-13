# ATTEND — bimodality measures for the TMB distribution reports
#
# Extracted from the former report 08 when it split into 07-tmb_distribution and
# 41_tcga_tmb_replication: both run the same three-criterion battery (dip test, Sarle's
# coefficient, BIC-selected mixture), 41 on the TCGA reference row, so keeping the
# definitions in one report would have meant copying them into the other.
#
# Hartigan's dip test itself lives in attend_classes.R::dip_bimodality(); this file adds
# the two measures that are not used anywhere else in the pipeline.
#
# No behaviour change: the bodies below are the former `bimod-helpers` chunk verbatim.

# --- Bimodality machinery used only in this report --------------------------------
# (Hartigan's dip test is reused from attend_classes.R::dip_bimodality(); the two
#  measures below are not reused elsewhere yet, so they live here.)

# Bias-corrected sample skewness / excess kurtosis (SAS/`moments` "type 2" convention,
# so the bimodality coefficient matches its published 5/9 benchmark). Computed inline
# to avoid adding a `moments` dependency.
.skew_kurt <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 4) return(c(skew = NA_real_, kurt = NA_real_))
  s <- sd(x); if (!is.finite(s) || s == 0) return(c(skew = NA_real_, kurt = NA_real_))
  z  <- (x - mean(x)) / s
  g1 <- (n / ((n - 1) * (n - 2))) * sum(z^3)                                  # skewness
  g2 <- ((n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))) * sum(z^4) -
        (3 * (n - 1)^2) / ((n - 2) * (n - 3))                                 # excess kurtosis
  c(skew = g1, kurt = g2)
}

# Sarle's bimodality coefficient b = (g1^2 + 1) / (g2 + 3(n-1)^2/((n-2)(n-3))).
# Benchmark: a uniform distribution gives b = 5/9 ≈ 0.555; b > 0.555 leans bimodal,
# b < 0.555 leans unimodal. (A heavy right tail alone inflates skewness and can push
# b up without a true second mode — hence we never rely on b alone.)
.bimodality_coef <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  sk <- .skew_kurt(x)
  if (any(!is.finite(sk))) return(NA_real_)
  (sk["skew"]^2 + 1) / (sk["kurt"] + (3 * (n - 1)^2) / ((n - 2) * (n - 3)))
}

# BIC-selected Gaussian mixture (1..3 components). Returns the chosen component count,
# the raw gap between the extreme component means, and that gap normalised by the
# vector's SD (the raw gap is NOT comparable across definitions, because the correction
# shifts the mut/Mb values; the SD-normalised gap is). Knit-safe without `mclust`.
.gmm_modes <- function(x) {
  x <- x[is.finite(x)]
  na <- list(G = NA_integer_, gap = NA_real_, gap_sd = NA_real_, mc = NULL)
  if (length(x) < 6 || !requireNamespace("mclust", quietly = TRUE)) return(na)
  mc <- tryCatch(mclust::Mclust(x, G = 1:3, verbose = FALSE), error = function(e) NULL)
  if (is.null(mc)) return(na)
  mns <- sort(as.numeric(mc$parameters$mean))
  gap <- if (length(mns) >= 2) mns[length(mns)] - mns[1] else 0
  s   <- sd(x)
  list(G = mc$G, gap = gap, gap_sd = if (is.finite(s) && s > 0) gap / s else NA_real_, mc = mc)
}

have_diptest <- requireNamespace("diptest", quietly = TRUE)
have_mclust  <- requireNamespace("mclust",  quietly = TRUE)
if (!have_diptest) message("`diptest` not installed — dip-test columns render as NA.")
if (!have_mclust)  message("`mclust` not installed — mixture columns render as NA.")

# --- shared derivations for report 07, 37 and 41 -----------------------------
#
# These were briefly written to disk by report 07 and read back by 37 and 41. That worked,
# but it imposed a BUILD ORDER on the reports for no reason: a shared derivation belongs in
# a function, not in an intermediate file. Both are pure — no I/O, no side effects — so the
# three reports can be built in any order and still get identical values.

#' One row per (patient, TMB definition), on the raw mut/Mb scale.
#'
#' Definitions come from tmb_defs_present(), so this cannot drift from the rest of the
#' pipeline. Definitions with fewer than `min_finite` finite values are dropped: the
#' ancestry-corrected column is all-NA until attend_ancestry$clinical_col is set, and a
#' distribution test needs values. The `definition` factor keeps correction order, so every
#' facet and axis reads normal -> ancestry corrected.
#'
#' @param master the joined master table
#' @param min_finite minimum finite values for a definition to be kept
#' @return tibble(pid, definition, tmb), with attr "finite_n" = finite count per definition
tmb_long <- function(master, min_finite = 10L) {
  def_cols <- tmb_defs_present(master)
  finite_n <- vapply(def_cols,
                     function(cc) sum(is.finite(suppressWarnings(as.numeric(master[[cc]])))),
                     integer(1))
  def_cols   <- def_cols[finite_n >= min_finite]
  def_labels <- tmb_def_labels(def_cols)

  out <- if (length(def_cols) == 0) {
    tibble::tibble(pid = character(0), definition = factor(), tmb = numeric(0))
  } else {
    master |>
      dplyr::select(pid, dplyr::all_of(unname(def_cols))) |>
      tidyr::pivot_longer(dplyr::all_of(unname(def_cols)),
                          names_to = "col", values_to = "tmb") |>
      dplyr::mutate(
        tmb        = suppressWarnings(as.numeric(tmb)),
        definition = factor(names(def_cols)[match(col, def_cols)],
                            levels = names(def_cols), labels = unname(def_labels))) |>
      dplyr::filter(is.finite(tmb)) |>
      dplyr::select(pid, definition, tmb)
  }
  attr(out, "finite_n") <- finite_n
  out
}

#' The three-criterion bimodality battery, one row per group.
#'
#' dip test (Hartigan), Sarle's bimodality coefficient, and a BIC-selected Gaussian
#' mixture. `verdict` is a majority vote: no single criterion is trustworthy on a
#' right-skewed distribution, which is the whole reason there are three.
#'
#' Report 41 calls this on the TCGA reference with `by = NULL` to get a single row it can
#' bind onto the ATTEND rows, which is why the grouping column is a parameter.
#'
#' @param d a data frame with a numeric `tmb` column
#' @param by name of the grouping column, or NULL for one row over all of `d`
tmb_battery <- function(d, by = "definition") {
  one <- function(x) {
    dp  <- dip_bimodality(x)
    bc  <- as.numeric(.bimodality_coef(x))
    gmm <- .gmm_modes(x)
    tibble::tibble(
      n           = length(x),
      dip_D       = dp$dip,
      dip_p       = dp$p,
      dip_bimodal = isTRUE(dp$bimodal),
      BC          = bc,
      BC_bimodal  = is.finite(bc) && bc > 0.555,
      gmm_G       = gmm$G,
      gmm_bimodal = is.finite(gmm$G) && gmm$G >= 2,
      gmm_gap_sd  = gmm$gap_sd)
  }
  res <- if (is.null(by)) one(d$tmb) else
    d |> dplyr::group_by(.data[[by]]) |>
      dplyr::group_modify(function(.x, ...) one(.x$tmb)) |> dplyr::ungroup()

  res |> dplyr::mutate(votes   = dip_bimodal + BC_bimodal + gmm_bimodal,
                       verdict = dplyr::if_else(votes >= 2, "bimodal", "unimodal"))
}
