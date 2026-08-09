# Plot helpers must return NULL rather than error on empty input — data/gistic/ is
# absent on any machine that has not run GISTIC, and the report must still knit.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  cat("ggplot2 absent — skipping (knit-safe by design)\n"); quit(status = 0)
}
source(file.path("code", "attend_classes.R"))
source(file.path("code", "attend_plots.R"))

peaks <- data.frame(peak_id = c("Amplification Peak 1", "Deletion Peak 1"),
                    descriptor = c("8q24.21", "10q23.31"),
                    direction = c("amp", "del"), chrom = c("8", "10"),
                    wide_start = c(1e6, 2e6), wide_end = c(2e6, 3e6),
                    q_value = c(1e-8, 1e-6), stringsAsFactors = FALSE)

freq <- data.frame(
  peak_id    = rep(peaks$peak_id, each = 2),
  scna_group = factor(rep(c("MMRd-high", "MMRp-high"), 2),
                      levels = attend_scna$group_levels),
  n = 9, n_altered = c(4, 5, 2, 3), freq = c(.44, .55, .22, .33),
  ci_lo = c(.14, .21, .03, .07), ci_hi = c(.79, .86, .60, .70),
  stringsAsFactors = FALSE)

p1 <- scna_mirror_plot(freq, peaks); stopifnot(inherits(p1, "ggplot"))
p2 <- scna_delta_plot(freq, peaks);  stopifnot(is.null(p2) || inherits(p2, "ggplot"))

sc <- c(A = 0.5, B = 0.2, C = 0.8)
g  <- factor(c("MMRd-high", "MMRd-high", "MMRp-high"),
             levels = attend_scna$group_levels)
p3 <- panel_score_plot(sc, g); stopifnot(inherits(p3, "ggplot"))

# Empty input -> NULL, not an error.
stopifnot(is.null(scna_mirror_plot(freq[0, ], peaks)))
stopifnot(is.null(panel_score_plot(numeric(0), factor(character(0)))))

cat("scna plot helpers: OK\n")
