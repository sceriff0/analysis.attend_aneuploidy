# Base-R unit test for add_response_class() / response_counts() (attend_classes.R §5b).
#
# The rule under test is the THIRD category. A landmark split on survival time has an
# obvious two-way reading -- past the cut is a responder, at-or-before is not -- and that
# reading is wrong for anyone whose follow-up merely STOPPED before the cut. Coding a
# patient censored at month 3 as a non-responder scores lost follow-up as treatment
# failure. This test pins the NA, because the bug it prevents is silent: the figures still
# render, the counts still add up, and only the responder arm is quietly wrong.
#
# Same library() shim as test_km_panel_counts.R, so attend_classes.R sources on a
# tidyverse-less machine (CLAUDE.md "Local dev caveat"). add_response_class() is pure base R.
local({
  real_library <- base::library
  assign("library", function(package, ...) {
    pkg <- as.character(substitute(package))
    if (length(find.package(pkg, quiet = TRUE)) == 0) return(invisible())
    real_library(pkg, character.only = TRUE, ...)
  }, envir = globalenv())
})
source(file.path("code", "attend_classes.R"))

cut <- attend_thresholds$response_months
stopifnot(is.numeric(cut), cut > 0)

#                     past cut   at cut    before cut  before cut  past cut   no time
#                     + event    + event   + event     censored    censored
df <- data.frame(
  gianlu__PFS_MONTHS = c(12,        cut,      3,          3,          24,        NA),
  gianlu__PFS_EVENT  = c("PD",      "PD",     "Death",    "No event", "No event", "PD"),
  stringsAsFactors = FALSE)

out <- add_response_class(df)
got <- as.character(out$response_class)

# 1. Past the landmark -> responder, EVEN WITH an event afterwards. The question is what
#    had happened by the landmark, not eventual outcome.
stopifnot(identical(got[1], "responder"))
stopifnot(identical(got[5], "responder"))

# 2. At or before the landmark WITH an event -> non-responder. The boundary is closed on
#    the non-responder side (time <= cut), so a patient progressing exactly at the cut is
#    not a responder.
stopifnot(identical(got[2], "non-responder"))
stopifnot(identical(got[3], "non-responder"))

# 3. THE POINT: censored before the landmark is unclassifiable, not a non-responder.
stopifnot(is.na(got[4]))

# 4. No survival time at all is unclassifiable too.
stopifnot(is.na(got[6]))

# Factor level order drives the x-axis and the fill palette in report 06.
stopifnot(identical(levels(out$response_class), c("responder", "non-responder")))
stopifnot(all(c("response_time", "response_event", "response_class") %in% names(out)))

# The cut is configuration, not a constant: raising it past month 12 must demote patient 1.
out20 <- add_response_class(df, thr = list(response_months = 20))
stopifnot(identical(as.character(out20$response_class)[1], "non-responder"),
          identical(as.character(out20$response_class)[5], "responder"))

# Missing survival columns give an all-NA class rather than an error, so a report still
# knits before attend_cols is pointed at the real column names (the have_cols() philosophy).
empty <- add_response_class(data.frame(pid = c("a", "b")))
stopifnot(nrow(empty) == 2, all(is.na(empty$response_class)))

# response_counts() must account for EVERY row: the two arms plus the two unclassifiable
# reasons partition the cohort, which is what makes it usable as a denominator.
if (length(find.package("tibble", quiet = TRUE)) > 0) {
  cnt <- response_counts(out)
  stopifnot(nrow(cnt) == 4, sum(cnt$n) == nrow(df))
  stopifnot(cnt$n[1] == 2, cnt$n[2] == 2, cnt$n[3] == 1, cnt$n[4] == 1)
} else {
  cat("test_response_class: tibble absent — response_counts() checks skipped\n")
}

cat("test_response_class: ALL PASS\n")
