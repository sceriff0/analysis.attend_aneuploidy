# Two denominator bugs in the recurrent-CNA figures (report 05):
#  (1) recurrent_arm_calls() folded ASCETS's LOWCOV (below breadth-of-coverage) into a
#      "neutral" catch-all, inflating the arm denominator with arms that were never
#      actually evaluated. ASCETS's own aneuploidy score, and this repo's
#      load_wes_results.R:485, both divide by sum(call != "LOWCOV") — the recurrence
#      barplot must match.
#  (2) segment_pileup() derived n_samples from the FILTERED segment table, so a sample
#      whose genome carried no surviving altered segment silently dropped out of the
#      denominator instead of counting as "profiled, not altered".
# These tests reproduce both bugs (task-1-brief.md) and lock in the fix.
source(file.path("code", "attend_classes.R"))

# ---- recurrent_arm_calls() -------------------------------------------------

# [1] LOWCOV excluded from the arm denominator: 4 samples, 2 AMP / 1 NEUTRAL / 1 LOWCOV
#     -> gain = 2/3 evaluable arms (0.667), not 2/4 (0.50).
calls1 <- tibble::tibble(ID = paste0("S", 1:4), `1q` = c("AMP", "AMP", "NEUTRAL", "LOWCOV"))
r1 <- recurrent_arm_calls(calls1)
stopifnot(abs(r1$freq$gain[r1$freq$arm == "1q"] - 2/3) < 1e-9)
stopifnot(abs(r1$freq$loss[r1$freq$arm == "1q"] - 0) < 1e-9)
stopifnot(abs(r1$freq$altered[r1$freq$arm == "1q"] - 2/3) < 1e-9)

# [2] NC, empty string, and real NA all behave as no-call (same as LOWCOV).
calls2 <- tibble::tibble(ID = paste0("S", 1:4), `2p` = c("AMP", "NC", "", NA))
r2 <- recurrent_arm_calls(calls2)
long2 <- r2$long[r2$long$arm == "2p", ]
stopifnot(is.na(long2$call[long2$ID == "S2"]))  # NC
stopifnot(is.na(long2$call[long2$ID == "S3"]))  # ""
stopifnot(is.na(long2$call[long2$ID == "S4"]))  # real NA
stopifnot(!is.na(long2$call[long2$ID == "S1"]) && long2$call[long2$ID == "S1"] == "gain")
stopifnot(r2$freq$n_evaluable[r2$freq$arm == "2p"] == 1)

# [3] An unrecognised token maps to NA and raises a single warning naming it.
calls3 <- tibble::tibble(ID = paste0("S", 1:3), `3p` = c("AMP", "WEIRDTOKEN", "NEUTRAL"))
warned <- FALSE
warn_msg <- NULL
r3 <- withCallingHandlers(
  recurrent_arm_calls(calls3),
  warning = function(w) { warned <<- TRUE; warn_msg <<- conditionMessage(w); invokeRestart("muffleWarning") }
)
stopifnot(warned)
stopifnot(grepl("WEIRDTOKEN", warn_msg, fixed = TRUE))
long3 <- r3$long[r3$long$arm == "3p", ]
stopifnot(is.na(long3$call[long3$ID == "S2"]))

# unrecognised tokens are collected ONCE per call, not once per token (only one warning).
calls3b <- tibble::tibble(ID = paste0("S", 1:4),
                          `3p` = c("FOO", "FOO", "BAR", "AMP"))
n_warnings <- 0
invisible(withCallingHandlers(
  recurrent_arm_calls(calls3b),
  warning = function(w) { n_warnings <<- n_warnings + 1; invokeRestart("muffleWarning") }
))
stopifnot(n_warnings == 1)

# [4] n_evaluable / n_lowcov are correct (from case [1]: 3 evaluable, 1 lowcov).
stopifnot(r1$freq$n_evaluable[r1$freq$arm == "1q"] == 3)
stopifnot(r1$freq$n_lowcov[r1$freq$arm == "1q"]    == 1)

# [5] min_evaluable_frac keeps a poorly-covered but highly-altered arm out of top_arms,
#     while it still appears in freq. "1p" is evaluable in only 2/10 samples (all AMP,
#     altered = 1.0 -- would rank #1 by frequency); "1q" is evaluable in 8/10, altered = 0.5.
ids5 <- paste0("S", 1:10)
calls5 <- tibble::tibble(
  ID  = ids5,
  `1p` = c("AMP", "AMP", rep("LOWCOV", 8)),
  `1q` = c(rep("AMP", 4), rep("NEUTRAL", 4), "LOWCOV", "LOWCOV")
)
r5 <- recurrent_arm_calls(calls5, n_top = 1)
stopifnot("1p" %in% r5$freq$arm)                 # still reported in freq
stopifnot(abs(r5$freq$altered[r5$freq$arm == "1p"] - 1.0) < 1e-9)
stopifnot(!("1p" %in% r5$top_arms))               # excluded: 2/10 < 0.5 min_evaluable_frac
stopifnot("1q" %in% r5$top_arms)                  # 8/10 >= 0.5, and now the top-ranked eligible arm

# [5b] min_evaluable_frac's denominator is DISTINCT ids, not nrow(calls) — load_ascets_calls()
#     applies the same id_strip collapse used elsewhere, so a resequenced sample can appear
#     under two rows. 4 DISTINCT samples (S1..S4); S1 is duplicated 3x as LOWCOV (4 rows for
#     S1), S2/S3/S4 each contribute one AMP row -> 7 total rows.
#     n_evaluable = 3 (S2, S3, S4). Row-count denominator (the bug): 3/7 = 0.429 < 0.5 ->
#     would wrongly exclude the arm. Distinct-id denominator (the fix): 3/4 = 0.75 >= 0.5 ->
#     correctly eligible, since only 4 real samples exist and 3 of them are evaluable.
calls5b <- tibble::tibble(
  ID  = c("S1", "S1", "S1", "S1", "S2", "S3", "S4"),
  `1p` = c("LOWCOV", "LOWCOV", "LOWCOV", "LOWCOV", "AMP", "AMP", "AMP")
)
r5b <- recurrent_arm_calls(calls5b, n_top = 1)
stopifnot(r5b$freq$n_evaluable[r5b$freq$arm == "1p"] == 3)
stopifnot(r5b$freq$eligible[r5b$freq$arm == "1p"])        # 3/4 distinct ids, not 3/7 rows
stopifnot("1p" %in% r5b$top_arms)

# ---- segment_pileup() -------------------------------------------------------

# [6] Pileup denominator counts samples PROFILED, not samples with a surviving altered
#     segment: 4 samples profiled, 2 carry a DUP segment -> gain = 0.5, not 1.0.
segs6 <- tibble::tibble(ID = c("S1", "S2"), Chromosome = "chr1",
                        Start = c(0, 0), End = c(5e5, 5e5), Type = "DUP", logRatio = 1.0)
arms6 <- tibble::tibble(chrom = "chr1", arm = c("1p", "1q"),
                        start = c(0, 1e6), end = c(1e6, 2e6))
p6 <- segment_pileup(segs6, arms6, bin = 1e6, min_abs = 0, n_samples = 4)
stopifnot(abs(p6$gain[1] - 0.5) < 1e-9)
stopifnot(attr(p6, "n_samples") == 4)

# [7] An explicit n_samples argument overrides the derived value.
#     segs6 has only 2 IDs (S1, S2), both with a DUP -> derived n_samples = 2, gain = 2/2 = 1.0.
#     Passing n_samples = 10 (e.g. the full cohort, most of it not in this segment table)
#     must instead give gain = 2/10 = 0.2.
p7_derived  <- segment_pileup(segs6, arms6, bin = 1e6, min_abs = 0)
p7_explicit <- segment_pileup(segs6, arms6, bin = 1e6, min_abs = 0, n_samples = 10)
stopifnot(abs(p7_derived$gain[1] - 1.0) < 1e-9)
stopifnot(attr(p7_derived, "n_samples") == 2)
stopifnot(abs(p7_explicit$gain[1] - 0.2) < 1e-9)
stopifnot(attr(p7_explicit, "n_samples") == 10)

# Denominator is also correctly derived from the UNFILTERED cnv_long when a sample's
# only segment is dropped by min_abs -- it must still count in the denominator.
segs7b <- tibble::tibble(ID = c("S1", "S2", "S3"), Chromosome = "chr1",
                         Start = c(0, 0, 0), End = c(5e5, 5e5, 5e5),
                         Type = "DUP", logRatio = c(1.0, 1.0, 0.01))
p7b <- segment_pileup(segs7b, arms6, bin = 1e6, min_abs = 0.5)
stopifnot(attr(p7b, "n_samples") == 3)            # S3's quiet/sub-threshold segment still counts
stopifnot(abs(p7b$gain[1] - 2/3) < 1e-9)          # only S1, S2 survive min_abs

# [8] NULL in -> NULL out holds for both functions (GC2).
stopifnot(is.null(recurrent_arm_calls(NULL)))
stopifnot(is.null(recurrent_arm_calls(tibble::tibble())))
stopifnot(is.null(segment_pileup(NULL, arms6)))
stopifnot(is.null(segment_pileup(segs6, NULL)))

# [9] Return shapes still carry every field named in GC4.
stopifnot(all(c("long", "freq", "top_arms", "wide_top") %in% names(r1)))
stopifnot(all(c("arm", "gain", "loss", "altered") %in% names(r1$freq)))
stopifnot(all(c("chrom", "bin_start", "gpos", "gain", "loss") %in% names(p6)))
stopifnot(!is.null(attr(p6, "chrom_offsets")))
stopifnot(!is.null(attr(p6, "n_samples")))

cat("All recurrent-CNA denominator properties verified.\n")
