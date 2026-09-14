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

## --- REGRESSION: the real matrix is UNSIGNED --------------------------------
## all_lesions.conf_*.txt stores a deletion call as +1/+2 (test_load_gistic_lesions.R).
## The old selection took "direction" from the sign of the group difference and tested
## `sign * call >= 1`, so a locus enriched in group B was never altered for anyone and
## the panel could only find group-A loci. Both deletion peaks below are unsigned: one
## enriched in A, one in B. Both must be selected AND both must move the score.
m <- rbind(matrix(c(1, 0), nrow = 4, ncol = 2, byrow = TRUE),
           matrix(c(0, 2), nrow = 4, ncol = 2, byrow = TRUE))
colnames(m) <- c("p_delA", "p_delB"); rownames(m) <- paste0("s", 1:8)
g  <- factor(rep(c("A", "B"), each = 4), levels = c("A", "B"))
pd <- c("del", "del")

sel <- select_panel_loci(m, g, n_select = 2, peak_dir = pd)
chk(identical(sort(colnames(m)[sel$idx]), c("p_delA", "p_delB")), "both peaks selected")
chk(identical(unname(sel$enrich[order(colnames(m)[sel$idx])]), c(1, -1)),
    "enrichment recorded per locus: p_delA -> A, p_delB -> B")
s <- .selected_panel_score(.peak_altered(m, pd), sel)
chk(all(s[1:4] == 1) && all(s[5:8] == 0),
    "a locus enriched in B contributes: A scores 1, B scores 0 (the old code gave A 0.5)")
res0 <- perm_test_selected_panel(m, g, n_select = 2, B = 200L, seed = 1, peak_dir = pd)
chk(abs(res0$stat - 1) < 1e-12, "statistic = mean |frequency difference| over selected loci")
chk(identical(sort(res0$enriched_in), c("A", "B")) && all(res0$peak_dir == "del"),
    "the output separates peak direction (del) from enrichment group (A / B)")

## --- a SIGNED matrix with peak directions: concordant events only ------------
sg <- rbind(matrix(c(2, -2), nrow = 4, ncol = 2, byrow = TRUE),
            matrix(c(0,  0), nrow = 4, ncol = 2, byrow = TRUE))
colnames(sg) <- c("p_amp", "p_del")
chk(all(.peak_altered(sg, c("amp", "del"), signed = TRUE)[1:4, ]), "signed: -2 at a del peak counts")
## No negative value left anywhere: signedness must come from the argument, not the data.
wrong <- sg; wrong[1:4, "p_del"] <- 2
chk(!any(.peak_altered(wrong, c("amp", "del"), signed = TRUE)[1:4, 2]),
    "signed: +2 at a del peak does NOT count (concordant events only)")
chk(all(.peak_altered(wrong, c("amp", "del"))[1:4, 2]),
    "unsigned default: +2 at a del peak IS a deletion call, as in all_lesions")

## --- a zero-difference peak is not selectable -------------------------------
flat <- cbind(m, p_flat = 1)
sel3 <- select_panel_loci(flat, g, n_select = 3)
chk(!("p_flat" %in% colnames(flat)[sel3$idx]), "a constant peak is excluded from selection")
chk(all(sel3$enrich != 0), "no selected locus carries enrichment 0")

## --- THE CALIBRATION CHECK: p must be uniform under the null ----------------
## Pure noise, random labels, nothing to find. 120 datasets, 40 samples, 60 peaks.
## The nested test must spread p over [0, 1]; the hoisted version must not.
naive_p <- function(mat, grp, n_select, B, seed) {
  # THE BUG, on purpose: selection ONCE, on the real labels, then permute the scores.
  A   <- .peak_altered(mat)
  sel <- .select_on_altered(A, grp, n_select)
  s   <- .selected_panel_score(A, sel)
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
chk(all(res$enriched_in[res$loci %in% paste0("pk", 1:10)] == "A"), "planted block enriched in A")
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
