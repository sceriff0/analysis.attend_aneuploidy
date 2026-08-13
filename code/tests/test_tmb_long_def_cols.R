# tmb_long() must expose the FILTERED key -> column map, aligned with the factor levels.
#
# This exists because of a real cluster failure: report 08's mode-validation chunk did
#   key_by_lab <- setNames(names(def_cols), unname(def_labels))
# where def_labels came from levels(long$definition) but def_cols was never assigned --
# it is a LOCAL inside tmb_long(), left behind when that function moved to
# attend_bimodality.R. The build died with "object 'def_cols' not found".
#
# The obvious repair, def_cols <- tmb_defs_present(master), is WRONG and would have failed
# silently instead of loudly: tmb_long() drops definitions under min_finite, so the
# unfiltered map is longer than the label vector the moment one definition is sparse --
# which is the normal state while tmb__TMB_nassar is all-NA. setNames() would then pair
# each label with the wrong column and every mode-validation table would silently describe
# the wrong TMB definition. Hence: expose the filtered map, and pin the alignment here.
#
# Same library() shim as test_km_panel_counts.R (no tidyverse meta-package locally), but
# dplyr/tidyr/tibble ARE present, so tmb_long() genuinely RUNS.
local({
  real_library <- base::library
  assign("library", function(package, ...) {
    pkg <- as.character(substitute(package))
    if (length(find.package(pkg, quiet = TRUE)) == 0) return(invisible())
    real_library(pkg, character.only = TRUE, ...)
  }, envir = globalenv())
})
source(file.path("code", "attend_classes.R"))
source(file.path("code", "attend_bimodality.R"))

if (length(find.package("dplyr", quiet = TRUE)) == 0 ||
    length(find.package("tidyr", quiet = TRUE)) == 0) {
  cat("test_tmb_long_def_cols: dplyr/tidyr absent — skipped\n"); quit(save = "no")
}

defs <- tmb_defs_present(data.frame(tmb__TMB_SCORE = 1, tmb__TMB_nassar = 1))
stopifnot(length(defs) == 2)   # guard: this test assumes both definitions are configured

n <- 20L
mk <- function(nassar) data.frame(pid = as.character(seq_len(n)),
                                  tmb__TMB_SCORE  = seq_len(n),
                                  tmb__TMB_nassar = nassar,
                                  stringsAsFactors = FALSE)

# --- both definitions populated -------------------------------------------------------
both <- tmb_long(mk(seq_len(n) + 0.5))
dc   <- attr(both, "def_cols")
stopifnot(!is.null(dc), length(dc) == 2)
# THE INVARIANT: one label per column, same order, so setNames() pairs them correctly.
stopifnot(identical(unname(tmb_def_labels(dc)), levels(both$definition)))
key_by_lab <- setNames(names(dc), levels(both$definition))
for (lab in levels(both$definition))
  stopifnot(identical(dc[[key_by_lab[[lab]]]], unname(dc)[match(lab, levels(both$definition))]))

# --- one definition all-NA (the CURRENT real state: nassar unset) ----------------------
# tmb_long() must drop it, and def_cols must shrink WITH the labels, not stay at 2.
one <- tmb_long(mk(NA_real_))
dc1 <- attr(one, "def_cols")
stopifnot(length(dc1) == 1, length(levels(one$definition)) == 1)
stopifnot(identical(unname(tmb_def_labels(dc1)), levels(one$definition)))
stopifnot(identical(unname(dc1), "tmb__TMB_SCORE"))
# The naive repair would have produced a length-2 map here — the silent-misalignment bug.
stopifnot(length(tmb_defs_present(mk(NA_real_))) != length(dc1))

# finite_n is deliberately PRE-filter (it reports what each definition had), so it stays
# at full length even when def_cols shrinks. Do not "fix" one to match the other.
stopifnot(length(attr(one, "finite_n")) == 2)

# --- nothing populated -----------------------------------------------------------------
none <- tmb_long(data.frame(pid = as.character(seq_len(n))))
stopifnot(length(attr(none, "def_cols")) == 0, length(levels(none$definition)) == 0)

cat("test_tmb_long_def_cols: ALL PASS\n")
