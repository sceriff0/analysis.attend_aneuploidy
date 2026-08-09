# Peak boundaries SHIFT between GISTIC runs (different cohorts -> different
# confidence intervals on peak location). Exact-coordinate matching would report
# every run's peaks as distinct and understate overlap, so union membership is
# defined by interval overlap, same chromosome, same direction.
source(file.path("code", "load_wes_results.R"))
source(file.path("code", "attend_classes.R"))

a <- data.frame(chrom = "8",  wide_start = 100, wide_end = 200, direction = "amp")
b <- data.frame(chrom = "8",  wide_start = 150, wide_end = 260, direction = "amp")
c_ <- data.frame(chrom = "8", wide_start = 300, wide_end = 400, direction = "amp")
d <- data.frame(chrom = "8",  wide_start = 150, wide_end = 260, direction = "del")
e <- data.frame(chrom = "10", wide_start = 150, wide_end = 260, direction = "amp")

stopifnot(peaks_overlap(a, b))            # overlapping, same chrom+direction
stopifnot(!peaks_overlap(a, c_))          # disjoint
stopifnot(!peaks_overlap(a, d))           # direction differs
stopifnot(!peaks_overlap(a, e))           # chromosome differs
stopifnot(peaks_overlap(a, data.frame(chrom = "8", wide_start = 200, wide_end = 300,
                                      direction = "amp")))  # 1 bp touch

# --- union across two synthetic run folders ---------------------------------
mk <- function(dirname, unique_name, descriptor, limits, s1) {
  dir.create(p <- file.path(tempdir(), dirname), showWarnings = FALSE)
  les <- data.frame(
    `Unique Name` = unique_name, Descriptor = descriptor,
    `Wide Peak Limits` = limits, `Peak Limits` = "na", `Region Limits` = "na",
    `q values` = 1e-6,
    `Residual q values after removing segments shared with higher peaks` = 1e-4,
    `Broad or Focal` = "Focal", `Amplitude Threshold` = "x",
    S1 = s1, check.names = FALSE)
  write.table(les, file.path(p, "all_lesions.conf_99.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  p
}
p1 <- mk("run_pooled", "Amplification Peak 1", "8q24", "chr8:100-200(probes 1:2)", 2)
p2 <- mk("run_mmrdhi", "Amplification Peak 1", "8q24", "chr8:150-260(probes 1:2)", 1)

cnv <- list(gistic = list(dir = "x", all_lesions_glob = "*all_lesions*.txt"),
            seg = list(id_strip = "-1TAD104|_tumor_only"))

u <- gistic_peak_union(c(pooled = p1, mmrd_high = p2), cnv)
stopifnot(nrow(u) == 2)                                  # both rows retained
stopifnot(identical(sort(unique(u$source)), c("mmrd_high", "pooled")))
stopifnot(length(unique(u$union_id)) == 1)               # they overlap -> one union peak
stopifnot("union_id" %in% names(u))

stopifnot(is.null(gistic_peak_union(c(a = file.path(tempdir(), "nope")), cnv)))

cat("peaks_overlap + gistic_peak_union: OK\n")
