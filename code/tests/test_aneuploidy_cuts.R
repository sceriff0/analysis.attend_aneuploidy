# add_aneuploidy_cuts(): the same patients, once per aneuploidy cut, with the class
# RE-DERIVED at each cut. The whole site draws both cuts from this one helper, so what
# needs pinning is that the second block is genuinely a different classification and not a
# relabelled copy of the first — a bug that would publish two identical panels under two
# cut labels and read as "the result does not depend on the threshold".
source(file.path("code", "attend_classes.R"))

fail <- character(0)
note <- function(...) fail <<- c(fail, paste0(...))

# Scores chosen to straddle BOTH cuts: 0.05 is low at either, 0.15 changes class between
# them, 0.25 is high at either, and NA must stay NA at both.
df <- data.frame(
  pid                    = c("a", "b", "c", "d"),
  aneu__aneuploidy_score = c(0.05, 0.15, 0.25, NA),
  gianlu__MMR_STATUS     = c("Deficient", "Intact", "Deficient", "Intact"),
  stringsAsFactors = FALSE
)

cuts <- aneu_cuts()
if (!identical(cuts, c(0.1, 0.2)))
  note("aneu_cuts() is ", paste(cuts, collapse = ", "), ", expected 0.1 and 0.2")

out <- add_aneuploidy_cuts(df, scna = TRUE)

## ---- shape: every patient appears once per cut -----------------------------------
if (nrow(out) != nrow(df) * length(cuts))
  note("add_aneuploidy_cuts() returned ", nrow(out), " rows, expected ",
       nrow(df) * length(cuts))
if (!all(c("as_cut", "as_cut_lab", "aneuploidy_class", "scna_group") %in% names(out)))
  note("missing output column(s): ",
       paste(setdiff(c("as_cut", "as_cut_lab", "aneuploidy_class", "scna_group"),
                     names(out)), collapse = ", "))
if (!identical(sort(unique(out$as_cut)), cuts))
  note("as_cut carries ", paste(sort(unique(out$as_cut)), collapse = ", "),
       ", expected ", paste(cuts, collapse = ", "))

## ---- the strip reads in cut order, not alphabetically ----------------------------
# Ordered so a facet row runs 0.1 then 0.2; a character column would sort the same way here
# but not once a third cut lands between them.
if (!is.factor(out$as_cut_lab) || !identical(levels(out$as_cut_lab),
                                             c("AS cut 0.1", "AS cut 0.2")))
  note("as_cut_lab levels are ",
       paste(levels(out$as_cut_lab), collapse = " | "), ", expected AS cut 0.1 | AS cut 0.2")

## ---- the class is RE-DERIVED, not copied -----------------------------------------
cls <- function(k) as.character(out$aneuploidy_class[out$as_cut == k])
if (!identical(cls(0.1), c("aneuploidy-low", "aneuploidy-high", "aneuploidy-high", NA)))
  note("at cut 0.1 the classes are ", paste(cls(0.1), collapse = ", "))
if (!identical(cls(0.2), c("aneuploidy-low", "aneuploidy-low", "aneuploidy-high", NA)))
  note("at cut 0.2 the classes are ", paste(cls(0.2), collapse = ", "))
# The patient at 0.15 is the whole point: high at 0.1, low at 0.2.
if (identical(cls(0.1), cls(0.2)))
  note("both cuts produced the SAME classification — the second block is a copy, so every ",
       "side-by-side panel would show one result twice under two labels")

## ---- the level order survives the re-derivation ----------------------------------
# Reference state first, as add_molecular_classes() guarantees; a facet must not flip the
# axis direction between the two cuts.
for (k in cuts) {
  lv <- levels(out$aneuploidy_class[out$as_cut == k])
  if (!identical(lv, c("aneuploidy-low", "aneuploidy-high")))
    note("at cut ", k, " aneuploidy_class levels are ", paste(lv, collapse = " | "))
}

## ---- scna_group follows the cut ---------------------------------------------------
# The MMR x AS grouping keys the Fig-1a column split and report 10's strata. Patient b is
# MMRp and moves high -> low, so its group must move with it.
grp <- function(k) as.character(out$scna_group[out$as_cut == k])
if (!identical(grp(0.1)[2], "MMRp-high") || !identical(grp(0.2)[2], "MMRp-low"))
  note("scna_group did not follow the cut: patient b is ", grp(0.1)[2], " at 0.1 and ",
       grp(0.2)[2], " at 0.2")
# NA score -> NA group at every cut, never a defaulted "low" (a default is not a measurement).
if (!all(is.na(grp(0.1)[4]), is.na(grp(0.2)[4])))
  note("a patient with no aneuploidy score was given a group")

## ---- the primary cut is unchanged by any of this ----------------------------------
# A bare add_molecular_classes() must still mean the PRIMARY cut: the master, the GISTIC
# run folders on disk and every unfaceted figure depend on it.
base_cls <- as.character(add_molecular_classes(df)$aneuploidy_class)
if (!identical(base_cls, cls(attend_thresholds$aneuploidy)))
  note("add_molecular_classes() no longer agrees with the primary cut ",
       attend_thresholds$aneuploidy)

## ---- collapsing back to one cut is one edit ---------------------------------------
one <- add_aneuploidy_cuts(df, cuts = attend_thresholds$aneuploidy)
if (nrow(one) != nrow(df)) note("a single-cut call returned ", nrow(one), " rows")

if (length(fail)) {
  cat("test_aneuploidy_cuts: ", length(fail), " failure(s)\n", sep = "")
  cat(paste0("  - ", fail, collapse = "\n"), "\n", sep = "")
} else {
  cat("test_aneuploidy_cuts: OK\n")
}
stopifnot(length(fail) == 0L)
