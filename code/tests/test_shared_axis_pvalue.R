# frac_display_y() — the shared axis, and the p-value that has to sit inside it.
#
# WHY THIS EXISTS. The IHC composition facets moved from scales = "free_y" to one shared
# axis per abundance tier. ggpubr derives stat_compare_means()'s default `label.y` from the
# PLOT-WIDE y range, not from the panel's own data, so on a shared axis every facet gets the
# same default, pinned to the global maximum. Measured on the four-facet repro below — one
# patient at ~10x the rest — all four p-values are computed correctly and drawn at y = 0.60
# while their boxes sit under 0.09: adrift at the ceiling, which on the page reads as no
# p-value at all, and the boxes occupy 13% of panel height. Two rounds of fixes went at the
# symptom ("restore the p-values") because the code LOOKED right: the labels were never
# missing, only misplaced, and nothing errors.
#
# So the check is a CALIBRATION, not a unit test. It asserts the property that failed —
# every label lands inside its panel, and the boxes fill the panel — against a dataset built
# to reproduce the failure. It fails against the regression (labels at the global max) and
# passes against frac_display_y().
#
# The numeric half is base R and always runs. The ggplot half needs ggplot2 + ggpubr and
# skips without them, so this file still runs in the bootstrap env.

source(file.path("code", "attend_plots.R"))

stopifnot(is.function(frac_display_y), is.function(shared_frac_y))

## --- the fixture: one extreme patient, everyone else in a narrow band -------
set.seed(1)
n  <- 12
d  <- do.call(rbind, lapply(c("CD8", "CD4", "Treg", "Bcell"), function(ct)
  data.frame(cell_type        = ct,
             aneuploidy_class = rep(c("low", "high"), each = n),
             .v               = c(rnorm(n, 0.05, 0.01), rnorm(n, 0.06, 0.01)),
             stringsAsFactors = FALSE)))
d$.v[1] <- 0.60                      # the patient that ranged the whole figure
d$.v[d$.v < 0] <- 0

y <- frac_display_y(d, ".v", c("aneuploidy_class", "cell_type"))

## The ceiling is the largest upper whisker, so it is nowhere near the outlier ...
stopifnot(y$ylim[2] < 0.2, y$ylim[2] > 0.07)
## ... the outlier is disclosed rather than dropped ...
## The wording is terse because it rides in a subtitle; what must hold is that a single
## off-scale point is COUNTED and singular, never silently dropped.
stopifnot(y$n_above == 1L, grepl("^; 1 point off scale", y$note))
## ... and the label sits inside the range, not on its edge.
stopifnot(y$label_y > 0, y$label_y < y$ylim[2])

## --- edge inputs must not produce a broken scale ----------------------------
## Flat at a non-zero level is NOT degenerate — it gets a real window, floored at 0.
flat <- data.frame(g = c("a", "a", "b", "b"), f = "x", v = c(1, 1, 1, 1))
yf <- frac_display_y(flat, "v", c("g", "f"))
stopifnot(!is.null(yf$ylim), yf$ylim[1] == 0, yf$ylim[2] > 1, yf$n_above == 0L)
## Nothing above the floor, and nothing finite at all, both fall back to ggplot's scaling.
zero  <- data.frame(g = c("a", "a", "b", "b"), f = "x", v = c(0, 0, 0, 0))
empty <- data.frame(g = character(0), f = character(0), v = numeric(0))
stopifnot(is.null(frac_display_y(zero,  "v", c("g", "f"))$ylim))
stopifnot(is.null(frac_display_y(empty, "v", c("g", "f"))$ylim))
## A NULL result must still be safe to pass straight into the axis helper.
stopifnot(length(shared_frac_y(frac_display_y(zero, "v", c("g", "f")))) == 1)
stopifnot(length(shared_frac_y(yf)) == 2)

## --- the calibration: labels inside the panel, boxes filling it -------------
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("ggpubr",  quietly = TRUE)) {

  # ggpubr must be ATTACHED, not just available: stat_compare_means() defers the label to
  # after_stat(create_p_label(...)), and that expression is resolved against the search
  # path, so a `ggpubr::` call alone dies with "could not find function create_p_label".
  library(ggpubr)

  p <- ggplot2::ggplot(d, ggplot2::aes(aneuploidy_class, .v)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::facet_wrap(~ cell_type) +
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format",
                               size = 3, label.y = y$label_y) +
    shared_frac_y(y)

  b   <- ggplot2::ggplot_build(p)
  lab <- b$data[[which(vapply(p$layers,
           function(l) inherits(l$stat, "StatCompareMeans"), logical(1)))]]
  rng <- b$layout$panel_params[[1]]$y.range

  # Every label is drawn, and every one is inside its panel.
  stopifnot(nrow(lab) == 4L, all(nzchar(as.character(lab$label))))
  stopifnot(all(lab$y > rng[1]), all(lab$y < rng[2]))

  # The boxes fill the panel instead of hugging the floor. Against the regression this
  # ratio is 0.13; the fix puts it above 0.6.
  stopifnot((0.09 - 0) / diff(rng) > 0.6)

  # coord_cartesian() zooms; it must NOT drop rows before stat_boxplot computes, or the
  # medians the figure reports are not the cohort's medians.
  bx <- b$data[[1]]
  stopifnot(nrow(bx) == 8L, all(is.finite(bx$middle)))
  full <- vapply(split(d$.v, paste(d$cell_type, d$aneuploidy_class)), stats::median, 0)
  stopifnot(max(abs(sort(bx$middle) - sort(unname(full)))) < 1e-9)
} else {
  message("test_shared_axis_pvalue.R: ggplot2/ggpubr absent — numeric half only.")
}

cat("test_shared_axis_pvalue.R: OK\n")
