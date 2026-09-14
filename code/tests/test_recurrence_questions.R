# Report 10's two questions, as tables — and the burden, concordance and run-size
# controls that make them readable.
#
# WHY THIS EXISTS. Four ways report 10 could answer its questions wrongly while every
# chunk rendered:
#  [1] Q1 read only the POOLED run's peaks, so a peak private to the nine-patient
#      MMRd-high run — the reason that run exists — never reached the Q1 table.
#  [2] Q1 called a peak "recurrent" when half the tumours carried it, a cut burden alone
#      clears. Recurrence is GISTIC's q on the stratum run.
#  [3] Q2's Spearman rho had no benchmark: shuffled splits of one population are the
#      ceiling "the same" can reach at n ~ 9.
#  [4] The panel permutation was described as holding burden fixed; it does not.
#
# Fixtures are UNSIGNED 0/1/2, like all_lesions.conf_*.txt: a deletion call is +1/+2.
# Base R only — sources attend_scna.R and lifts peaks_overlap() out of the loader
# without running the loader's library() calls.
source(file.path("code", "attend_scna.R"))
local({
  ex <- parse(file.path("code", "load_wes_results.R"))
  for (x in as.list(ex))
    if (is.call(x) && identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name("peaks_overlap")))
      eval(x, globalenv())
})
stopifnot(is.function(peaks_overlap))

fails <- 0L
chk <- function(ok, what) {
  cat(if (isTRUE(ok)) "  [PASS] " else "  [FAIL] ", what, "\n", sep = "")
  if (!isTRUE(ok)) fails <<- fails + 1L
}
lev <- c("MMRp-low", "MMRp-high", "MMRd-low", "MMRd-high")
pk <- function(id, desc, dir, chrom, s, e, q)
  data.frame(peak_id = id, descriptor = desc, direction = dir, chrom = chrom,
             wide_start = s, wide_end = e, q_value = q, stringsAsFactors = FALSE)

## --- [1] + [2] Q1 --------------------------------------------------------------
# Pooled run: 8q24 amp and 10q23 del. MMRd-high run: 8q24 amp (overlapping, different
# bounds) and a PRIVATE 19q12 amp the pooled run never called, plus a 13q14 del at q 0.4.
pooled_pk <- pk(c("Amplification Peak 1", "Deletion Peak 1"), c("8q24.21", "10q23.31"),
                c("amp", "del"), c("8", "10"), c(100, 500), c(200, 600), c(1e-6, 1e-3))
strat_pk  <- pk(c("Amplification Peak 1", "Amplification Peak 2", "Deletion Peak 1"),
                c("8q24.13", "19q12", "13q14.2"), c("amp", "amp", "del"),
                c("8", "19", "13"), c(150, 10, 40), c(250, 20, 50), c(0.01, 0.05, 0.4))
ids  <- c(paste0("d", 1:4), paste0("p", 1:4))
grp  <- factor(rep(c("MMRd-high", "MMRp-high"), each = 4), levels = lev)
pmat <- matrix(c(2, 0,  1, 1,  1, 0,  0, 2,     # d1-d4: 8q24 in 3/4
                 1, 1,  0, 1,  0, 1,  0, 0),    # p1-p4: 8q24 in 1/4; del 10q23 is +1 (unsigned)
               nrow = 8, byrow = TRUE, dimnames = list(ids, pooled_pk$peak_id))
smat <- matrix(c(1, 1, 0,  1, 1, 0,  0, 1, 1,  1, 1, 0), nrow = 4, byrow = TRUE,
               dimnames = list(ids[1:4], strat_pk$peak_id))
loo <- data.frame(peak_id = strat_pk$peak_id, retained_frac = c(1, 0.5, 0))

q1 <- q1_recurrence_table(list(peaks = pooled_pk, mat = pmat), list(peaks = strat_pk, mat = smat),
                          grp, fdr = 0.25, loo = loo)
chk("19q12" %in% q1$locus, "[1] a peak private to the stratum run reaches the Q1 table")
r19 <- q1[q1$locus == "19q12", ]
chk(r19$found_in == "MMRd-high run only" && is.na(r19$pooled_peak_id) &&
      r19$calls_from == "MMRd-high run" && r19$freq == 1,
    "[1] its frequency comes from the stratum run, and the row says so")
r8 <- q1[q1$locus == "8q24.13", ]
chk(r8$pooled_peak_id == "Amplification Peak 1" && r8$calls_from == "pooled run" &&
      r8$n == 4 && r8$n_altered == 3,
    "[1] an overlapping peak is scored from the POOLED run's calls, target patients only")
chk(!q1$recurrent[q1$locus == "13q14.2"], "[2] a stratum peak at q > fdr is not recurrent")
chk(!q1$recurrent[q1$locus == "10q23.31"] && q1$found_in[q1$locus == "10q23.31"] == "pooled run only",
    "[2] a pooled-only peak is listed but never called recurrent, however frequent")
chk(q1$loo_retained_frac[q1$locus == "19q12"] == 0.5, "leave-one-out retention is carried per locus")
chk(identical(q1$recurrent[1:2], c(TRUE, TRUE)), "recurrent loci sort first")

## --- Q2: conditional on Q1 -----------------------------------------------------
q2 <- q2_same_loci_table(q1, list(peaks = pooled_pk, mat = pmat), grp)
chk(setequal(q2$locus, c("8q24.13", "19q12")), "Q2 is restricted to Q1's recurrent loci")
chk(q2$freq_ref[q2$locus == "8q24.13"] == 0.25 && q2$diff[q2$locus == "8q24.13"] == 0.5,
    "Q2 scores the reference stratum on the pooled run")
chk(is.na(q2$freq_ref[q2$locus == "19q12"]) && grepl("not scorable", q2$note[q2$locus == "19q12"]),
    "a locus with no pooled peak is kept, NA, with the reason")

## --- [3] concordance null: identical populations reach the ceiling -------------
set.seed(5)
np <- 40L; n1 <- 9L; n2 <- 12L
prob <- runif(np, 0.05, 0.9)
same <- t(replicate(n1 + n2 + 10L, rbinom(np, 1, prob)))
colnames(same) <- paste0("pk", seq_len(np))
g3 <- factor(c(rep("MMRd-high", n1), rep("MMRp-high", n2), rep("MMRp-low", 10L)), levels = lev)
cn <- concordance_null(same, g3, B = 300)
chk(cn$null_lo < cn$rho && cn$null_hi > cn$null_lo, "[3] same-population rho falls inside the shuffled-split range")
flip <- same; flip[g3 == "MMRp-high", ] <- rbinom(n2 * np, 1, rep(1 - prob, each = n2))
cf <- concordance_null(flip, g3, B = 300)
chk(cf$rho < cf$null_lo, "[3] an inverted profile falls below the ceiling's 2.5% quantile")
chk(is.finite(cn$ref_rho) && cn$floor_group == "MMRp-low", "[3] the floor reference is reported")
chk(is.null(concordance_null(same[, 1:2], g3)), "fewer than 3 peaks -> NULL, not a rho of nothing")

## --- [4] burden control -------------------------------------------------------
# MMRp-high is altered at EVERY peak twice as often: a pure burden difference. The panel
# difference must sit inside the random-panel null, not outside it.
set.seed(9)
np <- 60L
bpk <- pk(paste0("pk", 1:np), paste0(rep(c("8q24", "17q12", "10q23", "9p21"), length.out = np), ".", 1:np),
          rep(c("amp", "amp", "del", "amp"), length.out = np), "1", 1:np, 1:np, 1e-3)
bpk$descriptor[1] <- "8q24.21"; bpk$descriptor[3] <- "10q23.31"
bg <- factor(rep(c("MMRp-high", "MMRd-high"), each = 30), levels = lev)
bm <- rbind(matrix(rbinom(30 * np, 1, 0.6), 30), matrix(rbinom(30 * np, 1, 0.3), 30))
colnames(bm) <- bpk$peak_id
panel2 <- data.frame(locus = c("MYC", "PTEN"), cytoband = c("8q24", "10q23"),
                     direction = c("amp", "del"), stringsAsFactors = FALSE)
pt <- perm_test_two_group(rowMeans(bm[, c("pk1", "pk3")]), bg, B = 500)
st <- panel_specificity_test(bm, bpk, bg, panel = panel2, B = 500)
chk(pt$p < 0.05, "[4] the label permutation DOES reject on burden alone (why the control exists)")
chk(st$p > 0.05 && st$n_amp == 1 && st$n_del == 1,
    sprintf("[4] the random-panel control does not (p = %.3f), matched on amp/del counts", st$p))
spec <- bm; spec[bg == "MMRd-high", c("pk1", "pk3")] <- 1L    # a real locus-specific effect
st2 <- panel_specificity_test(spec, bpk, bg, panel = panel2, B = 500)
chk(st2$p < 0.05, sprintf("[4] a locus-specific difference is detected (p = %.4f)", st2$p))

## --- run-size-aware overlap -----------------------------------------------------
u <- data.frame(source = c(rep("pooled", 10), "mmrd_high", "mmrd_high"),
                union_id = c(paste0("U", 1:10), "U1", "U2"), stringsAsFactors = FALSE)
ps <- peak_set_structure(u)
chk(ps$jaccard["pooled", "mmrd_high"] == 0.2 && ps$overlap["pooled", "mmrd_high"] == 1,
    "a small run that is a subset: Jaccard 0.2, overlap coefficient 1")
chk(identical(unname(ps$n_peaks[c("pooled", "mmrd_high")]), c(10L, 2L)), "peak counts per run reported")

if (fails > 0L) { cat("test_recurrence_questions: ", fails, " FAILED\n", sep = ""); quit(status = 1L) }
cat("test_recurrence_questions: ALL PASS\n")
