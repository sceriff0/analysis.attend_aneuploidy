# Base-R unit test; shim library() so attend_classes.R (library(tidyverse)) sources
# on a tidyverse-less machine. km_panel_counts is pure base R so it RUNS here.
local({
  real_library <- base::library
  assign("library", function(package, ...) {
    pkg <- as.character(substitute(package))
    if (length(find.package(pkg, quiet = TRUE)) == 0) return(invisible())
    real_library(pkg, character.only = TRUE, ...)
  }, envir = globalenv())
})
source(file.path("code", "attend_classes.R"))

df <- data.frame(
  gianlu__PFS_MONTHS = c(10, 20, 5, 30, 15, 8, 12, 40),
  gianlu__PFS_EVENT  = c("PD","No event","Death","No event","PD","PD","No event","Death"),
  aneuploidy_class   = c("aneuploidy-high","aneuploidy-high","aneuploidy-low","aneuploidy-low",
                         "aneuploidy-high","aneuploidy-high","aneuploidy-low","aneuploidy-low"),
  MMR_class          = c("Deficient","Deficient","Deficient","Deficient",
                         "Intact","Intact","Intact","Intact"),
  stringsAsFactors = FALSE)

res <- km_panel_counts(df, facets = c("aneuploidy_class","MMR_class"),
                       min_n = 3, min_events = 2)
stopifnot(is.data.frame(res), nrow(res) == 4,
          all(c("aneuploidy_class","MMR_class","n","n_events","small_n") %in% names(res)))
# Deficient × aneuploidy-high: rows 1,2 -> n=2, events PD+No event = 1
r1 <- res[res$aneuploidy_class=="aneuploidy-high" & res$MMR_class=="Deficient", ]
stopifnot(r1$n == 2, r1$n_events == 1, r1$small_n == TRUE)   # n=2 < min_n=3
# Intact × aneuploidy-high: rows 5,6 -> n=2, events PD+PD = 2
r2 <- res[res$aneuploidy_class=="aneuploidy-high" & res$MMR_class=="Intact", ]
stopifnot(r2$n == 2, r2$n_events == 2, r2$small_n == TRUE)
# NA rows dropped; missing facet column -> empty df
stopifnot(nrow(km_panel_counts(df, facets = "NOPE")) == 0)

cat("test_km_panel_counts: ALL PASS\n")
