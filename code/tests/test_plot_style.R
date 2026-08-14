# Figure-style enforcement for analysis/*.Rmd and the helper layer.
#
# Companion to test_rmd_style.R, and parse-only for the same reason: it reads the .Rmd as
# text and never evaluates a chunk, so it runs on a machine with no tidyverse and no cohort
# data (see CLAUDE.md "Local dev caveat") exactly as it does on the cluster. Base R only.
#
# It exists because the pipeline had drifted onto SIX palettes at once — Okabe-Ito,
# ColorBrewer Paired, ColorBrewer Set1/RdBu, matplotlib tab10, seaborn deep, and survminer's
# base-R colour names — and five themes. Nothing failed; every report knitted. The drift was
# only visible by reading ten files side by side, which is precisely the job a test should do.
#
# Six rules:
#   [1] no literal hex in a report — colour comes from the palette, not from the call site
#   [2] no raw ggplot theme_*() in a report — attend_theme() is the one theme
#   [3] no base-R colour NAME as a colour/fill constant ("red", "firebrick", "turquoise3")
#   [4] every geom_boxplot in a report also draws its points
#   [5] the highlight-group colours collide with no semantic palette colour
#   [6] km_facet()'s default palette is not a hardcoded pair of colour names

source(file.path("code", "attend_plots.R"))   # base-R-sourceable: no top-level library()

reports <- sort(list.files(file.path("analysis"), pattern = "\\.Rmd$", full.names = TRUE))
stopifnot(length(reports) > 0)

fail <- character(0)
note <- function(...) fail <<- c(fail, paste0(...))

# ---- helpers ---------------------------------------------------------------

# TRUE for every line inside a fenced block, so the prose rules never look at R code and
# the code rules never look at prose. Same shape as test_rmd_style.R::in_any_fence().
in_fence <- function(lines) {
  inside <- logical(length(lines)); open <- FALSE
  for (i in seq_along(lines)) {
    if (grepl("^```", lines[i])) { inside[i] <- TRUE; open <- !open; next }
    inside[i] <- open
  }
  inside
}

# Code lines only, with whole-line R comments dropped: a comment may legitimately NAME an
# old colour while explaining why it was replaced, and that must not trip rules [1] or [3].
code_lines <- function(lines) {
  keep <- in_fence(lines) & !grepl("^\\s*#", lines) & !grepl("^```", lines)
  data.frame(n = which(keep), txt = lines[keep], stringsAsFactors = FALSE)
}

# ---- [1] no literal hex in a report ----------------------------------------
for (f in reports) {
  cl <- code_lines(readLines(f, warn = FALSE))
  hit <- grepl('"#[0-9A-Fa-f]{6}"', cl$txt)
  if (any(hit))
    note(basename(f), ": literal hex in a report at line(s) ",
         paste(cl$n[hit], collapse = ", "),
         " — take the colour from attend_pal()/the semantic palettes in attend_plots.R, ",
         "so one edit moves every figure that uses it.")
}

# ---- [2] one theme ---------------------------------------------------------
# theme() (the modifier) is fine and expected; theme_<something>() is a competing base theme.
for (f in reports) {
  cl  <- code_lines(readLines(f, warn = FALSE))
  hit <- grepl("\\btheme_[a-z]+\\s*\\(", cl$txt) & !grepl("\\battend_theme\\s*\\(|\\btheme_pub\\s*\\(", cl$txt)
  if (any(hit))
    note(basename(f), ": raw ggplot theme at line(s) ", paste(cl$n[hit], collapse = ", "),
         " — reports use attend_theme(). Reports render to the workflowr HTML site, and a ",
         "second base theme changes type size and grid policy for that page alone.")
}

# ---- [3] no base-R colour names as constants -------------------------------
# Checked against grDevices::colors() rather than a hand-written denylist, so "present" or
# "mapped" (legend titles, factor levels) can never be mistaken for a colour. Greys, black,
# white and NA stay legal: they are the neutral/reference ink every figure needs.
r_colours <- grDevices::colors()
neutral   <- grepl("^(grey|gray)[0-9]*$|^(black|white|transparent)$", r_colours)
banned    <- r_colours[!neutral]
for (f in reports) {
  cl <- code_lines(readLines(f, warn = FALSE))
  for (i in seq_len(nrow(cl))) {
    m <- regmatches(cl$txt[i], gregexpr('(colour|color|fill)\\s*=\\s*"([^"]+)"', cl$txt[i]))[[1]]
    if (!length(m)) next
    vals <- sub('.*"([^"]+)"$', "\\1", m)
    bad  <- vals[vals %in% banned]
    if (length(bad))
      note(basename(f), ":", cl$n[i], ": base-R colour name(s) ",
           paste(sQuote(bad), collapse = ", "),
           " — these belong to no palette. Use attend_pal()/attend_neutral or a semantic ",
           "palette so the figure stays colour-blind-safe and consistent with its neighbours.")
  }
}

# ---- [4] a boxplot always shows its points ---------------------------------
# ATTEND is a ~40-patient cohort; split by MMR, aneuploidy class or responder status a group
# of four is routine. Every call site passes outlier.shape = NA (to avoid double-plotting),
# so a box with no point layer does not merely hide the outliers — it renders quartiles from
# a handful of values with nothing on the figure to say how few. attend_box() is the fix and
# carries its own points; a bare geom_boxplot() must supply them within a few lines.
POINT_LAYERS <- "geom_jitter|geom_point|highlight_points|aneu_point_layers|geom_quasirandom|geom_beeswarm"
for (f in reports) {
  lines <- readLines(f, warn = FALSE)
  for (i in grep("geom_boxplot\\s*\\(", lines)) {
    if (grepl("^\\s*#", lines[i])) next
    window <- lines[i:min(i + 6L, length(lines))]
    if (!any(grepl(POINT_LAYERS, window)))
      note(basename(f), ":", i, ": geom_boxplot() with no point layer within 6 lines — ",
           "prefer attend_box(), which draws the points and swaps the box for a median ",
           "crossbar wherever a group falls under n = ", attend_min_box_n, ".")
  }
}

# ---- [5] the highlight overlay must be visible on every box it lands on -----
# attend_classes.R is NOT base-R sourceable (it loads tidyverse), so the highlight colours are
# read as text. They are drawn ON TOP of boxes filled with the semantic palettes, so a highlight
# sharing a hex with any of them is invisible on exactly the box it exists to mark — which is
# what happened when polipo was #CC79A7 (= non-responder) and cohort was #E69F00 (= MMRd).
classes_src <- readLines(file.path("code", "attend_classes.R"), warn = FALSE)
hl_block <- grep("^attend_highlight\\s*<-", classes_src)
if (length(hl_block) == 1L) {
  tail_src <- classes_src[hl_block:min(hl_block + 40L, length(classes_src))]
  tail_src <- tail_src[!grepl("^\\s*#", tail_src)]
  hl_hex <- unlist(regmatches(tail_src, gregexpr('color\\s*=\\s*"#[0-9A-Fa-f]{6}"', tail_src)))
  hl_hex <- toupper(sub('.*"(#[0-9A-Fa-f]{6})"$', "\\1", hl_hex))
  semantic <- toupper(unique(c(attend_mmr_cols, attend_aneu_cols,
                               attend_tp53_cols, attend_resp_cols)))
  clash <- intersect(hl_hex, semantic)
  if (length(clash))
    note("attend_classes.R: highlight-group colour(s) ", paste(clash, collapse = ", "),
         " are also semantic palette colours — the highlight point disappears into any box ",
         "filled with the same hex. Use a hue no semantic palette claims.")
  if (!length(hl_hex))
    note("attend_classes.R: could not read any highlight-group colour — ",
         "rule [5] is silently passing; fix the parser or the attend_highlight block.")
} else {
  note("attend_classes.R: attend_highlight block not found — rule [5] cannot run.")
}

# ---- [6] km_facet() takes its palette from the project -----------------------
km_line <- grep("palette\\s*=", classes_src, value = TRUE)
km_line <- km_line[!grepl("^\\s*#", km_line)]
if (length(km_line)) {
  bad <- km_line[grepl(paste0('"(', paste(banned, collapse = "|"), ')"'), km_line)]
  if (length(bad))
    note("attend_classes.R: km_facet() default palette still names base-R colours (",
         trimws(bad[1]), ") — every survival curve in report 04 would be the one figure ",
         "family that cannot be read against the rest of the site.")
}

# ---- report ----------------------------------------------------------------
if (length(fail)) {
  cat("test_plot_style: FAIL\n"); cat(paste0("  - ", fail, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("test_plot_style: ALL PASS (", length(reports), " reports checked)\n", sep = "")
