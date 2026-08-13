# Task 2: genome ordering and amplitude tiers for the recurrent-CNA figures (report 06).
# Field convention (Beroukhim 2010; Zack 2013; GenVisR::cnFreq) is that recurrence is
# displayed ALONG THE GENOME, and that low-level gain/loss is separated from high-level
# amplification/homozygous deletion. Two new pieces:
#   order_arms_genomic() — reorders an arm-name vector into genomic order (a pure
#     transform; does not touch recurrent_arm_calls()'s own most-altered-first sort, GC4).
#   segment_pileup()'s gain_low/gain_high/loss_low/loss_high — split the existing gain/loss
#     columns (unchanged meaning, GC4) into amplitude tiers by |logRatio| >= high_abs.
source(file.path("code", "attend_classes.R"))

# ---- order_arms_genomic() --------------------------------------------------

# [1] Numeric chromosome order beats string order: "10p" sorts after "9q", before "11p"
#     (a plain string sort would put "10p" before "2p" entirely, and before "9q" too).
o1 <- order_arms_genomic(c("9q", "10p", "11p"))
stopifnot(identical(o1, c("9q", "10p", "11p")))
o1b <- order_arms_genomic(c("11p", "10p", "9q", "2p"))
stopifnot(identical(o1b, c("2p", "9q", "10p", "11p")))

# [2] Tolerates a "chr" prefix and mixed case; X sorts after 22, Y sorts after X.
o2 <- order_arms_genomic(c("chrXp", "22q", "CHR1p", "1q", "Yp", "chr1P"))
stopifnot(identical(o2, c("CHR1p", "chr1P", "1q", "22q", "chrXp", "Yp")))
# (two names both parse as "1p" — order_arms_genomic() is a pure reorder, so both survive,
#  and are tied on sort key; a stable sort keeps their relative input order, same as the
#  unparseable-name guarantee below.)

# [3] Unparseable names sort last, in original relative order, without erroring.
o3 <- order_arms_genomic(c("foo", "1p", "bar", "2q", "not_an_arm"))
stopifnot(identical(o3, c("1p", "2q", "foo", "bar", "not_an_arm")))

# order_arms_genomic() is a pure reorder: same length, same multiset of elements, no
# invented or dropped arms (ambiguity resolution #1).
in_arms <- c("2q", "1p", "chr17p", "Yq", "weird", "9q")
out_arms <- order_arms_genomic(in_arms)
stopifnot(length(out_arms) == length(in_arms))
stopifnot(setequal(out_arms, in_arms))
stopifnot(identical(sort(out_arms), sort(in_arms)))  # same multiset (all elements distinct here)

# ---- segment_pileup() amplitude tiers --------------------------------------

arms <- tibble::tibble(chrom = "chr1", arm = c("1p", "1q"),
                       start = c(0, 1e6), end = c(1e6, 2e6))

# [4] |logRatio| = 2 lands in gain_high, not gain_low (default high_abs = 1.0).
segs4 <- tibble::tibble(ID = "S1", Chromosome = "chr1", Start = 0, End = 5e5,
                        Type = "DUP", logRatio = 2)
p4 <- segment_pileup(segs4, arms, bin = 1e6, n_samples = 1)
stopifnot(abs(p4$gain_high[1] - 1) < 1e-9)
stopifnot(abs(p4$gain_low[1]  - 0) < 1e-9)

# [5] gain == gain_low + gain_high when no sample has both tiers in the same bin: one
#     low-magnitude sample, one high-magnitude sample, distinct ids.
segs5 <- tibble::tibble(ID = c("S1", "S2"), Chromosome = "chr1", Start = c(0, 0),
                        End = c(5e5, 5e5), Type = "DUP", logRatio = c(0.3, 2))
p5 <- segment_pileup(segs5, arms, bin = 1e6, n_samples = 2)
stopifnot(abs(p5$gain[1] - 1) < 1e-9)                        # 2/2 have a gain
stopifnot(abs(p5$gain_low[1]  - 0.5) < 1e-9)                 # S1 only
stopifnot(abs(p5$gain_high[1] - 0.5) < 1e-9)                 # S2 only
stopifnot(abs((p5$gain_low[1] + p5$gain_high[1]) - p5$gain[1]) < 1e-9)

# [6] A single sample with BOTH a low- and a high-magnitude gain segment in the same bin
#     counts ONCE, in gain_high — gain stays 1/n_samples (2 rows, 1 distinct id), not
#     2/n_samples, and gain_low + gain_high == gain still holds (ambiguity resolution #4).
segs6 <- tibble::tibble(ID = c("S1", "S1", "S2"), Chromosome = "chr1",
                        Start = c(0, 0, 0), End = c(5e5, 5e5, 5e5),
                        Type = "DUP", logRatio = c(0.3, 2, 0.1))
p6 <- segment_pileup(segs6, arms, bin = 1e6, n_samples = 2)
stopifnot(abs(p6$gain[1] - 1) < 1e-9)                         # 2/2 (S1, S2), not 3/2
stopifnot(abs(p6$gain_high[1] - 0.5) < 1e-9)                  # S1 only, once
stopifnot(abs(p6$gain_low[1]  - 0.5) < 1e-9)                  # S2 only
stopifnot(abs((p6$gain_low[1] + p6$gain_high[1]) - p6$gain[1]) < 1e-9)

# [7] NA logRatio counts as low tier.
segs7 <- tibble::tibble(ID = "S1", Chromosome = "chr1", Start = 0, End = 5e5,
                        Type = "DUP", logRatio = NA_real_)
p7 <- segment_pileup(segs7, arms, bin = 1e6, n_samples = 1)
stopifnot(abs(p7$gain_low[1]  - 1) < 1e-9)
stopifnot(abs(p7$gain_high[1] - 0) < 1e-9)

# [4b] Boundary: abs(logRatio) == high_abs lands in the HIGH tier (code uses >=, not >).
segs4b <- tibble::tibble(ID = "S1", Chromosome = "chr1", Start = 0, End = 5e5,
                         Type = "DUP", logRatio = 1.0)
p4b <- segment_pileup(segs4b, arms, bin = 1e6, n_samples = 1)
stopifnot(abs(p4b$gain_high[1] - 1) < 1e-9)

# same properties on the loss side, and with a custom high_abs threshold.
segs7b <- tibble::tibble(ID = c("S1", "S2"), Chromosome = "chr1", Start = c(0, 0),
                         End = c(5e5, 5e5), Type = "DEL", logRatio = c(-0.4, -1.5))
p7b <- segment_pileup(segs7b, arms, bin = 1e6, n_samples = 2, high_abs = 1)
stopifnot(abs(p7b$loss[1] - 1) < 1e-9)
stopifnot(abs(p7b$loss_low[1]  - 0.5) < 1e-9)
stopifnot(abs(p7b$loss_high[1] - 0.5) < 1e-9)
stopifnot(abs((p7b$loss_low[1] + p7b$loss_high[1]) - p7b$loss[1]) < 1e-9)

# [8] NULL in -> NULL out (GC2); GC4 fields (gain, loss, and the pre-existing shape) still
# present alongside the new tier columns.
stopifnot(is.null(segment_pileup(NULL, arms)))
stopifnot(is.null(segment_pileup(segs4, NULL)))
stopifnot(all(c("chrom", "bin_start", "gpos", "gain", "loss",
                "gain_low", "gain_high", "loss_low", "loss_high") %in% names(p6)))
stopifnot(!is.null(attr(p6, "chrom_offsets")))
stopifnot(!is.null(attr(p6, "n_samples")))

cat("All genome-order and amplitude-tier pileup properties verified.\n")
