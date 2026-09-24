# load_crispr_screens() — the Wu et al. (Immunity 2025) resister / sensitizer catalogue.
#
# WHY THIS EXISTS. data/crispr_t_cell_screens_wu2025.tsv is a derived reference table that
# nothing re-reads by eye. Four numbers in it are stated in the paper, so they are the only
# available check that the conversion still lines up with the source:
#
#   519 resisters and 877 sensitizers at >= 2 independent measurements   (Results, p. 2865)
#   108 resisters and 90 sensitizers overlapping their own recurrent set (Figures 1D, 1E)
#
# code/fetch_crispr_screens.R already stop()s on 519/877. It does NOT check 108/90, and that
# is the gap this closes: the workbook holds those two flags in ADJACENT column blocks
# (column M is "...DELETED genes" for resisters, column T is "...AMPLIFIED genes" for
# sensitizers), so a column swap leaves 519/877 untouched while inverting the direction of
# every claim the report makes. Same shape as the .peak_altered() sign bug.
#
# Also pinned: the evidence cutoff is a real ARGUMENT, not a relabelling of a pre-filtered
# file; both-list genes are dropped and counted; and the four gene records that sit on the
# workbook's second header row survive (a header-inference reader loses them).
#
# Skips cleanly when the TSV is absent, so it runs before the data lands.

fails <- 0L
chk <- function(ok, what) {
  cat(if (isTRUE(ok)) "  [PASS] " else "  [FAIL] ", what, "\n", sep = "")
  if (!isTRUE(ok)) fails <<- fails + 1L
}

tsv <- file.path("data", "crispr_t_cell_screens_wu2025.tsv")
if (!file.exists(tsv)) {
  cat("test_crispr_catalogue: SKIP — ", tsv, " not present\n", sep = ""); quit(status = 0L)
}

# Base R read, deliberately: this test must not depend on the loader's own parsing choices,
# or it would verify the loader against itself.
d <- utils::read.delim(tsv, stringsAsFactors = FALSE, colClasses = "character")
d$n_measurements <- suppressWarnings(as.integer(d$n_measurements))

chk(all(c("gene", "role", "n_measurements", "wu2025_dp_overlap") %in% names(d)),
    "the TSV carries gene, role, n_measurements, wu2025_dp_overlap")
chk(setequal(unique(d$role), c("resister", "sensitizer")),
    "role is exactly {resister, sensitizer}")

## --- the four published numbers ---------------------------------------------
at2 <- function(role) sum(d$role == role & d$n_measurements >= 2L, na.rm = TRUE)
chk(at2("resister")   == 519L, paste0("519 resisters at >= 2 measurements (got ", at2("resister"), ")"))
chk(at2("sensitizer") == 877L, paste0("877 sensitizers at >= 2 measurements (got ", at2("sensitizer"), ")"))

ov <- function(role) sum(d$role == role & d$n_measurements >= 2L & d$wu2025_dp_overlap == "YES",
                         na.rm = TRUE)
# THE CHECK fetch_crispr_screens.R's own stop() does not make. 108 is Figure 1D (deleted x
# resister); 90 is Figure 1E (amplified x sensitizer). If these two ever swap, 519/877 still
# pass and every directional statement in report 13 silently inverts.
chk(ov("resister")   == 108L, paste0("108 resisters carry Wu's DP-deleted overlap flag (got ",
                                     ov("resister"), ")"))
chk(ov("sensitizer") ==  90L, paste0("90 sensitizers carry Wu's DP-amplified overlap flag (got ",
                                     ov("sensitizer"), ")"))
chk(ov("resister") > ov("sensitizer"),
    "the resister overlap is the larger of the two — orientation sanity, independent of value")

## --- both-list genes exist, and are exactly the documented count ------------
both <- intersect(d$gene[d$role == "resister"   & d$n_measurements >= 2L],
                  d$gene[d$role == "sensitizer" & d$n_measurements >= 2L])
chk(length(both) == 105L, paste0("105 symbols are in BOTH lists at >= 2 (got ", length(both), ")"))
# Wu et al. never address both-list genes. Theirs are disjoint only on the OVERLAP flags,
# which is what makes their own directional claim safe — and what ATTEND must not assume.
both_flagged <- intersect(d$gene[d$role == "resister"   & d$wu2025_dp_overlap == "YES"],
                          d$gene[d$role == "sensitizer" & d$wu2025_dp_overlap == "YES"])
chk(length(both_flagged) == 0L,
    "no gene carries Wu's overlap flag in both directions (their 108 and 90 are disjoint)")

## --- the workbook's second-header-row genes survived the conversion ---------
## Rows 57 and 58 of the sheet are the second study block's title and header in columns A-G
## while simultaneously holding ordinary gene data in the gene blocks. A reader that infers
## a header for the whole grid drops these four records.
boundary <- c("KCTD5", "MLST8", "PIGS", "PIGU")
chk(all(boundary %in% d$gene),
    paste0("the four second-header-row genes survive (", paste(boundary, collapse = ", "), ")"))

## --- no duplicate symbol within a role --------------------------------------
dup <- any(duplicated(d[, c("gene", "role")]))
chk(!dup, "no gene symbol is duplicated within a role")

## --- the cutoff is an ARGUMENT, not a pre-filtered file ---------------------
## The unfiltered blocks must still be present, or `min_measurements` would be decorative:
## a file already cut at >= 2 would return the same answer for every value.
chk(sum(d$role == "resister")   == 2828L,
    paste0("the resister block is unfiltered in the file (2,828 rows, got ",
           sum(d$role == "resister"), ")"))
chk(sum(d$role == "sensitizer") == 2961L,
    paste0("the sensitizer block is unfiltered in the file (2,961 rows, got ",
           sum(d$role == "sensitizer"), ")"))
chk(at2("resister") < sum(d$role == "resister"),
    "raising min_measurements actually removes rows")

## --- the loader honours its arguments ---------------------------------------
if (requireNamespace("readr", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE) &&
    requireNamespace("here", quietly = TRUE)) {
  suppressPackageStartupMessages(source(file.path("code", "attend_classes.R")))
  suppressPackageStartupMessages(source(file.path("code", "load_crispr.R")))
  c2 <- load_crispr_screens(min_measurements = 2L, drop_both_list = TRUE)
  c1 <- load_crispr_screens(min_measurements = 1L, drop_both_list = TRUE)
  chk(nrow(c1) > nrow(c2), "min_measurements = 1 returns more genes than 2")
  chk(sum(c2$role == "resister") == 519L - length(both),
      "drop_both_list removes the 105 from the kept set, not from the count of 519")
  chk(length(attr(c2, "both_list")) == 105L,
      "the dropped both-list genes are returned as an attribute, so they can be printed")
  kept <- load_crispr_screens(min_measurements = 2L, drop_both_list = FALSE)
  chk(sum(kept$role == "resister") == 519L,
      "with drop_both_list = FALSE the raw 519 is recoverable")
  chk("wu_overlap" %in% names(c2) && !"wu2025_dp_overlap" %in% names(c2),
      "the ambiguous source column is renamed at load time")
} else {
  cat("  [SKIP] loader checks need readr/dplyr/here\n")
}

if (fails > 0L) {
  cat("test_crispr_catalogue: ", fails, " FAILED\n", sep = ""); quit(status = 1L)
}
cat("test_crispr_catalogue: ALL PASS\n")
