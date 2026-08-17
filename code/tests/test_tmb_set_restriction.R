# Base-R unit test. attend_classes.R calls library(tidyverse) at its top; on a
# machine without tidyverse (dev bootstrap, R 4.6) that blocks source(). The TMB
# config + tmb_defs_present() are pure base R, so install a library() shim that
# no-ops ONLY for packages that are not installed (and delegates to the real
# library() when they ARE — so this file behaves identically on the cluster).
local({
  real_library <- base::library
  assign("library", function(package, ...) {
    pkg <- as.character(substitute(package))
    if (length(find.package(pkg, quiet = TRUE)) == 0) return(invisible())
    # delegate with the character form: base::library() uses NSE, so passing the
    # promise `package` would make its own substitute() see the symbol `package`,
    # not the real name. Use the already-resolved `pkg` string + character.only.
    real_library(pkg, character.only = TRUE, ...)
  }, envir = globalenv())
})
source(file.path("code", "attend_classes.R"))  # sourced under the shim

## robust/matched are intermediates (computed, never plotted); the plotted set is
## the ancestry-corrected column only.
stopifnot(identical(unname(attend_cols$tmb_intermediate),
                    c("tmb__TMB_robust","tmb__TMB_matched")))
stopifnot(identical(unname(attend_cols$tmb_set), "tmb__TMB_nassar"),
          names(attend_cols$tmb_set) == "nassar")
## tmb_set_full is GONE — nothing may reintroduce a third plotted definition.
stopifnot(is.null(attend_cols$tmb_set_full))

## a master carrying every TMB column
df <- data.frame(pid = 1:2,
                 tmb__TMB_SCORE   = c(5, 12),
                 tmb__TMB_robust  = c(4, 11),
                 tmb__TMB_matched = c(4, 10),
                 tmb__TMB_nassar  = c(3, 9))

## THE only result: normal + ancestry corrected, in that order. There is no
## `full` argument any more, so no report can widen this.
d <- tmb_defs_present(df)
stopifnot(identical(unname(d), c("tmb__TMB_SCORE","tmb__TMB_nassar")),
          names(d)[1] == "upstream", names(d)[2] == "nassar")
stopifnot(!"full" %in% names(formals(tmb_defs_present)))

## even on a master carrying robust/matched, they are never returned
stopifnot(!any(c("tmb__TMB_robust","tmb__TMB_matched") %in% tmb_defs_present(df)))

## reader-facing labels: "nassar" must never reach the page
stopifnot(identical(tmb_def_labels(d),
                    c("TMB (normal)", "TMB (ancestry corrected)")))
## unmapped definitions fall back to the column name rather than erroring
stopifnot(identical(tmb_def_labels(c(mystery = "tmb__TMB_X")), "tmb__TMB_X"))

## track names: syntactically safe, and the internal name never appears in either form
stopifnot(identical(tmb_def_track_names(d), c("TMB_normal", "TMB_ancestry_corrected")),
          identical(tmb_def_track_names(d), make.names(tmb_def_track_names(d))))
stopifnot(!any(grepl("nassar|upstream", c(tmb_def_labels(d), tmb_def_track_names(d)),
                     ignore.case = TRUE)))

## the gene panels: MSH3 in the MMR union, Wnt genes deliberately OUT of it
stopifnot(all(c("PMS2","MLH1","MSH2","MSH6","MSH3") %in% attend_gene_panel),
          !any(c("CTNNB1","APC") %in% attend_gene_panel),
          identical(attend_gene_panel_wnt, c("CTNNB1","APC")),
          identical(attend_genes_per_gene, c(attend_gene_panel, attend_gene_panel_wnt)))

## highlight groups: `polipo` and nothing else. The retired `cohort` group marked patients
## holding IHC/imaging data in yellow -- data availability, not biology, on ~16 of ~40 patients,
## so it read as a second competing fill rather than a mark. Availability lives on the master
## as the in_* columns. test_figure_system.R carries the full reasoning and the colour pin.
stopifnot(identical(names(attend_highlight$groups), "polipo"))

cat("test_tmb_set_restriction: ALL PASS\n")
