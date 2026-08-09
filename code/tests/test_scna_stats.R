# Three properties that matter more than the point estimates:
#  (1) Wilson score CIs [U5], because 4/9 must display as an interval (~17-69%),
#      not as "44%"; Wilson is preferred over Clopper-Pearson at small n;
#  (2) permutation p uses the (1+k)/(B+1) correction, so p is never exactly 0 —
#      an uncorrected 0/B would be reported as p=0, which is not a claim the data
#      can support at n=9;
#  (3) Families A and B are corrected SEPARATELY — pooling would spend the
#      confirmatory panel's power on genome-wide discovery multiplicity.
source(file.path("code", "attend_classes.R"))
source(file.path("code", "attend_scna.R"))

mat <- matrix(c(2, 0,
                1, 0,
                0, 1,
                0, 0),
              nrow = 4, byrow = TRUE,
              dimnames = list(paste0("S", 1:4), c("P1", "P2")))
grp <- factor(c("MMRd-high", "MMRd-high", "MMRp-high", "MMRp-high"),
              levels = attend_scna$group_levels)

ft <- scna_freq_table(mat, grp)
r  <- ft[ft$peak_id == "P1" & ft$scna_group == "MMRd-high", ]
stopifnot(r$n == 2, r$n_altered == 2, r$freq == 1)
stopifnot(r$ci_lo < 1, r$ci_hi <= 1, r$ci_lo > 0)   # Wilson interval, not a point estimate
stopifnot(all(c("ci_lo", "ci_hi") %in% names(ft)))

set.seed(1)
x <- c(rep(0.1, 9), rep(0.6, 20))
g <- factor(c(rep("MMRd-high", 9), rep("MMRp-high", 20)))
pt <- perm_test_two_group(x, g, B = 2000)
stopifnot(pt$p > 0)                            # (1+k)/(B+1) — never exactly zero
stopifnot(pt$p < 0.01)
stopifnot(abs(pt$stat - (0.1 - 0.6)) < 1e-9)

p <- c(0.01, 0.02, 0.03)
stopifnot(all(family_adjust(p, "A") == p.adjust(p, "holm")))
stopifnot(all(family_adjust(p, "B") == p.adjust(p, "BH")))
stopifnot(!identical(family_adjust(p, "A"), family_adjust(p, "B")))

# [U1] the FDR/Holm denominator is the FULL union/panel count, including NA
# (untested) entries — NOT stats::p.adjust()'s own default, which silently drops
# NA before sizing the correction (n = 2, the smaller/wrong denominator here).
# family_adjust() must count the NA row toward n (n = 3).
p_na <- c(0.001, 0.01, NA)
stopifnot(identical(round(unname(family_adjust(p_na, "B")), 6),
                    round(unname(p.adjust(p_na, "BH", n = 3)), 6)))
stopifnot(!identical(round(unname(family_adjust(p_na, "B")), 6),
                     round(unname(p.adjust(p_na, "BH", n = 2)), 6)))

cat("scna_freq_table + perm_test_two_group + family_adjust: OK\n")
