# add_molecular_classes() returns ORDERED FACTORS with the reference state FIRST, and
# as_label() maps every DATA spelling in the pipeline onto ONE DISPLAY spelling.
#
# Both exist because of the same audit (2026-09-14). Before it:
#
#  * The classes were plain character vectors, so every axis, legend and facet strip got
#    ggplot's ALPHABETICAL order — which puts "aneuploidy-high" before "aneuploidy-low",
#    i.e. the ALTERED end on the left, against the direction attend_aneu_cols is built for
#    (light blue quiet -> red altered). Report 04 then re-factored locally to high-first
#    while reports 06 and 11 built their own low-first, so one variable read left-to-right
#    in opposite directions on different pages of one site.
#  * The same class was spelled three ways for aneuploidy ("aneuploidy-high", "aneu-high",
#    "AS High") and four ways for TP53 ("TP53-abnormal", "mut", "mutant", and the hyphenless
#    display form) — two of the TP53 pair in ADJACENT figures in report 07. The COLOURS were
#    right throughout, because attend_aneu_cols and attend_tp53_cols key every spelling.
#    That permissive key set is precisely what kept the drift invisible, so the guard has to
#    check the LABELS, not the hexes.
#
# DATA spellings must NOT change: add_scna_group() derives its groups with
# grepl("high", ...) and scna_group_token() turns those into GISTIC folder names on disk.
# The round-trip assertion below is what pins that.

source(file.path("code", "attend_classes.R"))
source(file.path("code", "attend_plots.R"))

fail <- character(0)
note <- function(...) fail <<- c(fail, paste0(...))

# ---- [1] reference-first ordered factors -----------------------------------
df <- data.frame(
  `aneu__aneuploidy_score` = c(0.05, 0.50, NA),
  `tmb__TMB_SCORE`         = c(1, 50, NA),
  `hrd__HRD_Score`         = c(1, 80, NA),
  `gianlu__MSI_STATUS`     = c("Stable (MSS)", "Instable (MSI)", NA),
  `gianlu__MMR_STATUS`     = c("Intact", "Deficient", NA),
  check.names = FALSE, stringsAsFactors = FALSE)
out <- add_molecular_classes(df)

want <- list(
  aneuploidy_class = c("aneuploidy-low", "aneuploidy-high"),
  TMB_class        = c("TMB-low", "TMB-high"),
  HRD_class        = c("HRD-low", "HRD-high"),
  MSI_class        = c("MSS", "MSI-high"),
  MMR_class        = c("Intact", "Deficient"))
for (nm in names(want)) {
  if (!is.factor(out[[nm]])) {
    note(nm, " is not a factor — ggplot would order it alphabetically, which puts the ",
         "ALTERED state first for every one of these.")
    next
  }
  if (!identical(levels(out[[nm]]), want[[nm]]))
    note(nm, " levels are ", paste(levels(out[[nm]]), collapse = " -> "),
         "; expected ", paste(want[[nm]], collapse = " -> "),
         " (reference state first, so the reader travels toward the altered end).")
}

# ---- [2] the DATA spellings survive the factor conversion ------------------
if (!identical(as.character(out$aneuploidy_class),
               c("aneuploidy-low", "aneuploidy-high", NA_character_)))
  note("aneuploidy_class no longer round-trips through as.character() — add_scna_group()'s ",
       "grepl(\"high\", ...) and scna_group_token()'s GISTIC folder names both read it.")

g <- add_scna_group(out)
if (!identical(as.character(g$scna_group), c("MMRp-low", "MMRd-high", NA_character_)))
  note("add_scna_group() no longer derives the same groups from the factor: got ",
       paste(as.character(g$scna_group), collapse = ", "))

# ---- [3] one display spelling per class ------------------------------------
disp <- list(
  # every DATA spelling that reaches a figure -> the one thing a reader should see
  `aneuploidy-high` = "AS High",   `aneuploidy-low` = "AS Low",
  `aneu-high`       = "AS High",   `aneu-low`       = "AS Low",
  `TP53-abnormal`   = "TP53 abnormal", `TP53-normal`  = "TP53 normal",
  `mut`             = "TP53 abnormal", `wt`           = "TP53 normal",
  `mutant`          = "TP53 abnormal", `wild-type`    = "TP53 normal",
  `Deficient`       = "MMRd",      `Intact`          = "MMRp",
  `MMR deficient`   = "MMRd",      `MMR proficient`  = "MMRp",
  `TMB-high`        = "TMB high",  `TMB-low`         = "TMB low")
for (k in names(disp)) {
  got <- as_label(k)
  if (!identical(got, disp[[k]]))
    note("as_label(\"", k, "\") is \"", got, "\", expected \"", disp[[k]],
         "\" — the same class must read the same way in a boxplot, a facet strip and a ",
         "heatmap annotation bar.")
}
# The stratum-prefixed form must keep its prefix.
if (!identical(as_label("MMRd aneuploidy-high"), "MMRd AS High"))
  note("as_label() dropped the stratum prefix on \"MMRd aneuploidy-high\".")
# Anything with no display spelling passes through untouched, so as_label() is safe as a
# blanket scale labeller on an axis that mixes classes.
if (!identical(as_label(c("something else", NA)), c("something else", NA_character_)))
  note("as_label() altered a value with no display spelling, or lost an NA.")

# ---- [4] the heatmap bar titles come from the shared constants -------------
if (!identical(unname(.fig1a_ann_titles[["Aneuploidy_hl"]]), attend_as_legend))
  note("the Fig-1a binary aneuploidy bar no longer takes its title from attend_as_legend — ",
       "it published the data-frame column name \"Aneuploidy_hl\" beside a boxplot legend ",
       "reading \"AS class\" for the same variable.")
for (k in c("MMR", "TP53"))
  if (!nzchar(.fig1a_ann_titles[[k]]))
    note("Fig-1a annotation ", k, " has no title constant.")

# ---- [4b] the palette keys the RAW clinical spelling too --------------------
# attend_mmr_cols hardcodes "Deficient"/"Intact" because attend_plots.R is sourced in a
# base-R bootstrap env and its order against attend_classes.R varies by report. This is the
# assertion that keeps the hardcoded pair in step with the configured token: without the raw
# keys, a fill scale over MMR_class matches nothing and every box renders in ggplot's NA grey
# with the legend silently dropped.
if (!attend_levels$mmr_deficient %in% names(attend_mmr_cols))
  note("attend_mmr_cols has no key for attend_levels$mmr_deficient (\"",
       attend_levels$mmr_deficient, "\") — a fill scale over MMR_class would render every ",
       "box in ggplot's NA grey and drop the legend.")
if (!identical(unname(attend_mmr_cols[[attend_levels$mmr_deficient]]),
               unname(attend_mmr_cols[["MMRd"]])))
  note("the raw clinical MMR key and \"MMRd\" resolve to different colours.")

# ---- [5] relabelling keeps each value paired with its colour ---------------
a <- data.frame(
  MMR           = factor(c("Intact", "Deficient"), levels = c("Intact", "Deficient")),
  TP53          = factor(c("wt", "mut"), levels = c("wt", "mut")),
  Aneuploidy_hl = factor(c("aneu-low", "aneu-high"), levels = c("aneu-low", "aneu-high")),
  stringsAsFactors = FALSE)
rl <- .fig1a_relabel_display(a, .fig1a_covariate_cols(a))
pairs <- list(c("TP53", "TP53 abnormal", unname(attend_tp53_cols[["mut"]])),
              c("TP53", "TP53 normal",   unname(attend_tp53_cols[["wt"]])),
              c("Aneuploidy_hl", "AS High", unname(attend_aneu_cols[["aneu-high"]])),
              c("Aneuploidy_hl", "AS Low",  unname(attend_aneu_cols[["aneu-low"]])))
for (pp in pairs) {
  # Defensive indexing: a broken relabel makes the key MISSING rather than wrong, and a
  # test that dies with "subscript out of bounds" on the very bug it guards reports nothing.
  cv  <- rl$cols[[pp[1]]]
  got <- if (!is.null(cv) && pp[2] %in% names(cv)) unname(cv[[pp[2]]]) else NA_character_
  if (!identical(got, pp[3]))
    note("after the display relabel, ", pp[1], " key \"", pp[2], "\" maps to ", got,
         " not ", pp[3], " — the values and the colour vector fell out of step, which is ",
         "the one way this relabel can silently mis-colour a bar.")
}

if (length(fail)) {
  cat("test_class_level_order: FAIL\n"); cat(paste0("  - ", fail, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("test_class_level_order: ALL PASS\n")
