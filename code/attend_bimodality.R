# ATTEND — bimodality measures for the TMB distribution reports
#
# Extracted from the former report 08 when it split into 36_tmb_distribution and
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
