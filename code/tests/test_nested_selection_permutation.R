# perm_test_selected_panel() — the data-driven (Family B) panel score.
#
# WHY THIS EXISTS. A data-driven panel is a legitimate analysis and an INVALID test if
# the selection sits outside the permutation. Selecting the most group-discriminating
# peaks uses the labels, so permuting the labels afterwards cannot reproduce the
# optimism the selection introduced: the null comes out far too narrow and p is near
# zero on data with no signal at all (Kriegeskorte et al., Nat Neurosci 2009 —
# "double dipping"). The fix is one line's POSITION: select_panel_loci() is called
# inside replicate(), on the permuted labels.
#
# That makes an ordinary unit test nearly useless here — both versions return a
# number in [0, 1] and the valid one is not "more correct" on any single dataset.
# What separates them is CALIBRATION over many null datasets, so that is what this
# checks: p must be roughly uniform when there is nothing to find. Hoisting the
# selection out of the replicate (the obvious "optimisation") fails this immediately.
#
# Base R only, so this runs in the bootstrap env without tidyverse.
source(file.path("code", "attend_scna.R"))

stopifnot(is.function(select_panel_loci), is.function(perm_test_selected_panel),
          is.function(perm_test_global_burden))

fails <- 0L
chk <- function(ok, what) {
  cat(if (isTRUE(ok)) "  [PASS] " else "  [FAIL] ", what, "\n", sep = "")
  if (!isTRUE(ok)) fails <<- fails + 1L
}

## --- direction is respected, not just magnitude -----------------------------
## Two peaks. Peak 1: group A amplified (+2), group B flat. Peak 2: group A deleted
## (-2), group B flat. A magnitude-only rule would score both the same; a directional
## rule must call peak 1 "amp" and peak 2 "del".
m <- rbind(matrix(c( 2, -2), nrow = 4, ncol = 2, byrow = TRUE),
           matrix(c( 0,  0), nrow = 4, ncol = 2, byrow = TRUE))
colnames(m) <- c("p_amp", "p_del"); rownames(m) <- paste0("s", 1:8)
g <- factor(rep(c("A", "B"), each = 4), levels = c("A", "B"))

sel <- select_panel_loci(m, g, n_select = 2)
chk(identical(sort(colnames(m)[sel$idx]), c("p_amp", "p_del")), "both peaks selected")
chk(identical(unname(sel$dir[order(colnames(m)[sel$idx])]), c(1, -1)),
    "direction taken from the sign: p_amp -> amp, p_del -> del")

## the -2 call counts as altered at the del-direction locus, and a +2 there would not
alt <- .directional_altered(m, sel$idx, sel$dir)
chk(all(alt[1:4, ]), "group A altered at both loci in their own directions")
chk(!any(alt[5:8, ]), "flat group A altered at neither")
flipped <- m; flipped[1:4, "p_del"] <- 2      # wrong-direction event
chk(!any(.directional_altered(flipped, sel$idx, sel$dir)[1:4, 2]),
    "a +2 at a del-direction locus does NOT count (concordant events only)")

## --- a zero-difference peak is not selectable -------------------------------
## sign(0) == 0 would make .directional_altered() test `0 >= 1` for every sample and
## contribute a constant column, silently shrinking the score's range.
flat <- cbind(m, p_flat = 0)
sel3 <- select_panel_loci(flat, g, n_select = 3)
chk(!("p_flat" %in% colnames(flat)[sel3$idx]), "a constant peak is excluded from selection")
chk(all(sel3$dir != 0), "no selected locus carries direction 0")

## --- THE CALIBRATION CHECK: p must be uniform under the null ----------------
## Pure noise, random labels, nothing to find. 120 datasets, 40 samples, 60 peaks.
## The nested test must spread p over [0, 1]; the hoisted version must not.
naive_p <- function(mat, grp, n_select, B, seed) {
  # THE BUG, on purpose: selection ONCE, on the real labels, then permute the scores.
  sel <- select_panel_loci(mat, grp, n_select)
  s   <- rowMeans(.directional_altered(mat, sel$idx, sel$dir), na.rm = TRUE)
  lv  <- levels(grp)
  obs <- mean(s[grp == lv[1]]) - mean(s[grp == lv[2]])
  set.seed(seed)
  null <- replicate(B, { gp <- sample(grp); mean(s[gp == lv[1]]) - mean(s[gp == lv[2]]) })
  (1 + sum(null >= obs)) / (B + 1)
}

set.seed(7)
R <- 120L; B <- 200L; n <- 40L; p <- 60L
p_nested <- numeric(R); p_naive <- numeric(R)
for (r in seq_len(R)) {
  mm <- matrix(sample(c(-2L, -1L, 0L, 0L, 1L, 2L), n * p, replace = TRUE), nrow = n)
  colnames(mm) <- paste0("pk", seq_len(p)); rownames(mm) <- paste0("s", seq_len(n))
  gg <- factor(rep(c("A", "B"), each = n / 2), levels = c("A", "B"))
  p_nested[r] <- perm_test_selected_panel(mm, gg, n_select = 8, B = B, seed = r)$p
  p_naive[r]  <- naive_p(mm, gg, n_select = 8, B = B, seed = r)
}

cat(sprintf("  null calibration over %d datasets:\n", R))
cat(sprintf("    nested : mean p = %.3f   P(p <= 0.05) = %.3f   min = %.4f\n",
            mean(p_nested), mean(p_nested <= 0.05), min(p_nested)))
cat(sprintf("    hoisted: mean p = %.3f   P(p <= 0.05) = %.3f   min = %.4f\n",
            mean(p_naive), mean(p_naive <= 0.05), min(p_naive)))

chk(mean(p_nested) > 0.30 && mean(p_nested) < 0.70,
    sprintf("nested mean p is central (%.3f in 0.30-0.70) -> roughly uniform", mean(p_nested)))
chk(mean(p_nested <= 0.05) < 0.20,
    sprintf("nested type-I error is not grossly inflated (%.3f)", mean(p_nested <= 0.05)))
chk(mean(p_naive) < 0.10,
    sprintf("hoisted version IS anti-conservative (mean p = %.3f) -> the bug is detectable",
            mean(p_naive)))
chk(mean(p_naive <= 0.05) > 3 * max(mean(p_nested <= 0.05), 0.01),
    "hoisted rejects far more often than nested on pure noise")

## --- and it still has power when the signal is real -------------------------
## Group A carries a genuine 10-peak amplified block; B does not.
set.seed(11)
mm <- matrix(sample(c(-1L, 0L, 0L, 1L), n * p, replace = TRUE), nrow = n)
colnames(mm) <- paste0("pk", seq_len(p)); rownames(mm) <- paste0("s", seq_len(n))
gg <- factor(rep(c("A", "B"), each = n / 2), levels = c("A", "B"))
mm[gg == "A", 1:10] <- 2L
res <- perm_test_selected_panel(mm, gg, n_select = 10, B = 400L, seed = 3)
cat(sprintf("  planted 10-peak block: p = %.4f, %d/%d planted peaks recovered\n",
            res$p, sum(res$loci %in% paste0("pk", 1:10)), 10L))
chk(res$p < 0.05, sprintf("nested test detects a real block (p = %.4f)", res$p))
chk(sum(res$loci %in% paste0("pk", 1:10)) >= 8, "selection recovers the planted peaks")
chk(all(res$dir[res$loci %in% paste0("pk", 1:10)] == "amp"), "planted block called amp")
chk(length(res$loo_frac) == length(res$loci) && all(res$loo_frac >= 0 & res$loo_frac <= 1),
    "leave-one-out selection stability reported, one fraction per selected locus")
chk(mean(res$loo_frac) > 0.8, sprintf("a real block is LOO-stable (mean %.2f)", mean(res$loo_frac)))

## --- the p floor is honest --------------------------------------------------
chk(res$p >= 1 / (400L + 1), "p cannot fall below the 1/(B+1) resolution floor")

## --- the global no-selection test still works -------------------------------
gb <- perm_test_global_burden(mm, gg, B = 400L, seed = 3)
chk(is.finite(gb$p) && gb$p >= 1 / (gb$B + 1), "perm_test_global_burden() returns a bounded p")
chk(gb$p < 0.05, sprintf("global burden also sees the planted block (p = %.4f)", gb$p))

if (fails > 0L) {
  cat("test_nested_selection_permutation: ", fails, " FAILED\n", sep = ""); quit(status = 1L)
}
cat("test_nested_selection_permutation: ALL PASS\n")
