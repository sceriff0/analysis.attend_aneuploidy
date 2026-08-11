# Task 6: segment_pileup()'s corrected denominator (Task 1) is opt-in via an explicit
# n_samples argument, and nothing in production passes one. The gap: a sample that
# contributes ZERO rows to cnv_long (a copy-number-quiet tumour whose ClassifyCNV file has
# no calls) cannot be counted by inspecting cnv_long's rows, no matter how early the count
# is taken relative to min_abs / direction filtering. load_cnv_data() is the only place
# that ever sees such a sample (one file per sample, including zero-row files), so it now
# attaches attr(out, "n_profiled") <- the count of DISTINCT derived sample ids among those
# files (6a, code/load_wes_results.R; not a raw file count — see .classifycnv_file_id() and
# test_classifycnv_file_id.R for why two files can de-duplicate to one sample id). This
# file tests the consumer side (6b): segment_pileup() prefers that attribute over the
# distinct-id count derived from cnv_long, with a floor guard and an explicit-argument
# override, recording which source won in attr(., "n_samples_source").
#
# 6a itself (load_cnv_data()) is NOT exercised here — it needs data.table, which is not
# installed on this machine. See the Task 6 report for how that was verified by reading
# the surrounding code paths instead.
source(file.path("code", "attend_classes.R"))

arms <- tibble::tibble(chrom = "chr1", arm = c("1p", "1q"),
                       start = c(0, 1e6), end = c(1e6, 2e6))

# Base fixture: 2 distinct ids present (S1, S2), both carrying a DUP segment.
base_segs <- tibble::tibble(ID = c("S1", "S2"), Chromosome = "chr1",
                            Start = c(0, 0), End = c(5e5, 5e5),
                            Type = "DUP", logRatio = 1.0)

# [1] attr(x, "n_profiled") <- 4 with only 2 ids present -> denominator is 4, not 2.
segs1 <- base_segs
attr(segs1, "n_profiled") <- 4
p1 <- segment_pileup(segs1, arms, bin = 1e6, min_abs = 0)
stopifnot(abs(p1$gain[1] - 2/4) < 1e-9)
stopifnot(attr(p1, "n_samples") == 4)
stopifnot(identical(attr(p1, "n_samples_source"), "profiled"))

# [2] No attribute present -> falls back to the observed distinct-id count (2).
p2 <- segment_pileup(base_segs, arms, bin = 1e6, min_abs = 0)
stopifnot(abs(p2$gain[1] - 2/2) < 1e-9)
stopifnot(attr(p2, "n_samples") == 2)
stopifnot(identical(attr(p2, "n_samples_source"), "observed"))

# [3] An explicit n_samples argument overrides a present attribute.
segs3 <- base_segs
attr(segs3, "n_profiled") <- 4
p3 <- segment_pileup(segs3, arms, bin = 1e6, min_abs = 0, n_samples = 10)
stopifnot(abs(p3$gain[1] - 2/10) < 1e-9)
stopifnot(attr(p3, "n_samples") == 10)
stopifnot(identical(attr(p3, "n_samples_source"), "argument"))

# [4] Floor guard: an attribute of 1 with 2 ids present does not shrink the denominator —
#     the observed count (2) wins, and the source is "observed" (the attribute was rejected,
#     not merely overridden).
segs4 <- base_segs
attr(segs4, "n_profiled") <- 1
p4 <- segment_pileup(segs4, arms, bin = 1e6, min_abs = 0)
stopifnot(abs(p4$gain[1] - 2/2) < 1e-9)
stopifnot(attr(p4, "n_samples") == 2)
stopifnot(identical(attr(p4, "n_samples_source"), "observed"))

# [5] NA attribute is ignored, falling back to observed.
segs5 <- base_segs
attr(segs5, "n_profiled") <- NA_integer_
p5 <- segment_pileup(segs5, arms, bin = 1e6, min_abs = 0)
stopifnot(abs(p5$gain[1] - 2/2) < 1e-9)
stopifnot(attr(p5, "n_samples") == 2)
stopifnot(identical(attr(p5, "n_samples_source"), "observed"))

# [6] NULL in -> NULL out still holds (GC2), attribute or not.
stopifnot(is.null(segment_pileup(NULL, arms)))
stopifnot(is.null(segment_pileup(base_segs, NULL)))
segs6 <- tibble::tibble(ID = character(0), Chromosome = character(0),
                        Start = numeric(0), End = numeric(0),
                        Type = character(0), logRatio = numeric(0))
stopifnot(is.null(segment_pileup(segs6, arms)))

cat("All pileup profiled-denominator properties verified.\n")
