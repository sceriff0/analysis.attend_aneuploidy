# The composite panel score is report 15's PRIMARY endpoint. Two failure modes it
# must not have:
#  (1) direction-blindness — a deletion at MYC scoring as serous-like;
#  (2) a fixed denominator — loci with no matching peak must leave the denominator,
#      not count as "not altered", which would deflate every sample equally and
#      silently shrink the effect.
source(file.path("code", "attend_classes.R"))
source(file.path("code", "attend_scna.R"))

peaks <- data.frame(
  peak_id    = c("Amplification Peak 1", "Deletion Peak 1", "Deletion Peak 2"),
  descriptor = c("8q24.21", "10q23.31", "8q24.21"),
  direction  = c("amp", "del", "del"),
  chrom      = c("8", "10", "8"),
  wide_start = c(1, 1, 1), wide_end = c(2, 2, 2),
  q_value    = c(1e-8, 1e-6, 1e-6),
  stringsAsFactors = FALSE
)
mat <- matrix(c(2, 0, 0,
                0, 2, 2,
                0, 0, 0),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("A", "B", "C"), peaks$peak_id))

panel <- data.frame(locus = c("MYC", "PTEN", "CCNE1"),
                    cytoband = c("8q24", "10q23", "19q12"),
                    direction = c("amp", "del", "amp"),
                    stringsAsFactors = FALSE)

mp <- match_panel_peaks(peaks, panel)
stopifnot(mp$peak_id[mp$locus == "MYC"]  == "Amplification Peak 1")  # amp, not the 8q24 del
stopifnot(mp$peak_id[mp$locus == "PTEN"] == "Deletion Peak 1")
stopifnot(is.na(mp$peak_id[mp$locus == "CCNE1"]))                    # 19q12 has no peak

s <- panel_score(mat, peaks, panel)
stopifnot(attr(s, "n_loci") == 2)          # CCNE1 leaves the denominator
stopifnot(s[["A"]] == 0.5)                 # MYC amp only
stopifnot(s[["B"]] == 0.5)                 # PTEN del only (8q24 DEL must not count as MYC)
stopifnot(s[["C"]] == 0)
stopifnot(identical(names(s), c("A", "B", "C")))

cat("match_panel_peaks + panel_score: OK\n")
