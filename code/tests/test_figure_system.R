# Base-R unit tests (no tidyverse) — runnable in the bootstrap env.
source(file.path("code", "attend_plots.R"))  # base-R-sourceable: no top-level library()

## --- semantic palettes exist and carry the promoted Fig-1a hexes ---
## The palette is Paul Tol's "vibrant" scheme, NOT Okabe-Ito: Okabe-Ito has no true red, so
## the aneuploidy scale could not be the light-blue -> red pair the project describes in prose
## without lying about the hex. See the palette block in attend_plots.R.
## The length pin is about attend_pal() and the k = 2..8 cluster bar, which draw from this
## vector POSITIONALLY and must get 8 distinct colours without recycling. So what has to hold
## is that the FIRST EIGHT slots are the cluster set and are unchanged; the aneuploidy pair
## appended after them (pale_blue/salmon, the restored ColorBrewer look) cannot resequence a
## cluster colour. A bare `== 8L` would have failed that append while pinning nothing extra.
stopifnot(is.character(ATTEND_PALETTE),
          length(ATTEND_PALETTE) >= 8L,
          identical(names(ATTEND_PALETTE)[1:8],
                    c("blue", "cyan", "teal", "orange", "red", "magenta", "grey", "black")),
          !any(duplicated(ATTEND_PALETTE)),
          !any(duplicated(ATTEND_PALETTE)),
          identical(unname(ATTEND_PALETTE[["grey"]]), "#999999"))

## Tol's own grey is #BBBBBB, which is indistinguishable from attend_neutral (#BFBFBF) and
## would recreate the bug this palette pass fixed — a non-semantic box reading as a semantic
## reference state. The deviation is deliberate; this pins it so a "fix the palette to match
## the paper" edit cannot land silently.
stopifnot(!"#BBBBBB" %in% ATTEND_PALETTE)

stopifnot(is.character(attend_mmr_cols),
          attend_mmr_cols[["MMRp"]] == "#999999",
          attend_mmr_cols[["MMRd"]] == "#EE7733",
          attend_mmr_cols[["MMR proficient"]] == "#999999",
          attend_mmr_cols[["MMR deficient"]] == "#EE7733")

stopifnot(identical(attend_aneu_low, "#D6E5F0"),   # pale blue = aneuploidy-low
          identical(attend_aneu_high, "#E3918F"),  # salmon    = aneuploidy-high
          attend_aneu_high != attend_tp53_cols[["mut"]],  # must not clash with the adjacent
          attend_aneu_high != attend_mmr_cols[["MMRd"]],  # Fig-1a annotation bars
          is.character(attend_aneu_cols),
          attend_aneu_cols[["aneuploidy-low"]]  == attend_aneu_low,
          attend_aneu_cols[["aneuploidy-high"]] == attend_aneu_high,
          attend_aneu_cols[["aneu-low"]]  == attend_aneu_low,
          attend_aneu_cols[["aneu-high"]] == attend_aneu_high,
          attend_aneu_cols[["MMRd aneuploidy-low"]]  == attend_aneu_low,
          attend_aneu_cols[["MMRd aneuploidy-high"]] == attend_aneu_high)

stopifnot(is.character(attend_tp53_cols),
          attend_tp53_cols[["TP53-normal"]]   == "#999999",
          attend_tp53_cols[["TP53-abnormal"]] == "#EE3377")

## --- the three rules the palette assignment rests on (attend_plots.R) -------------
## [A] grey is the reference state, everywhere. If a "normal"/"proficient"/wild-type level
##     ever picks up a hue, the reader loses the one cue that says "do not look here".
REFERENCE_GREY <- "#999999"
stopifnot(attend_mmr_cols[["MMRp"]]        == REFERENCE_GREY,
          attend_tp53_cols[["TP53-normal"]] == REFERENCE_GREY,
          attend_tp53_cols[["wt"]]          == REFERENCE_GREY)

## [B] concepts that can share a figure must not share a hue. The Fig-1a annotation stack
##     puts MMR, aneuploidy and TP53 side by side, so their ALTERED states are mutually
##     distinct. Response reuses the TP53 hue on purpose (report 06 is the only report that
##     colours by response and it has no TP53) — if that ever stops being true, this is the
##     test that has to change, not the figure.
altered <- c(attend_mmr_cols[["MMRd"]], attend_aneu_high, attend_tp53_cols[["mut"]])
stopifnot(length(unique(altered)) == 3L,
          !attend_aneu_low %in% altered)

## [C] `blue` (#0077BB) is RESERVED for the highlight overlay and no semantic palette may
##     claim it. Highlight points are drawn ON TOP of boxes filled with these colours, so a
##     highlight sharing a hex is invisible on exactly the box it exists to mark. Keeping one
##     hue unclaimed is what makes that guarantee cheap — and it is why the retired `cohort`
##     group had to reach for yellow, the least visible overlay in the set.
HIGHLIGHT_BLUE <- "#0077BB"
stopifnot(identical(unname(ATTEND_PALETTE[["blue"]]), HIGHLIGHT_BLUE),
          !HIGHLIGHT_BLUE %in% c(attend_mmr_cols, attend_aneu_cols,
                                 attend_tp53_cols, attend_resp_cols))

## Response must not reuse EITHER aneuploidy pole: reports 05 and 06 draw the same
## cell-type panels keyed on different variables, so a shared colour would make two
## different groupings look like one.
stopifnot(!attend_resp_cols[["responder"]]     %in% c(attend_aneu_low, attend_aneu_high),
          !attend_resp_cols[["non-responder"]] %in% c(attend_aneu_low, attend_aneu_high),
          attend_resp_cols[["responder"]] != attend_resp_cols[["non-responder"]],
          attend_resp_cols[["responder"]] == "#009988")

## The neutral single-group fill must not carry ANY semantic meaning — this is the bug that
## motivated the palette pass: #A6CEE3 was simultaneously "aneuploidy-low" and "no grouping".
semantic <- unique(c(attend_mmr_cols, attend_aneu_cols, attend_tp53_cols, attend_resp_cols))
stopifnot(is.character(attend_neutral), !attend_neutral %in% semantic)

## Every semantic colour is drawn from the one palette (attend_neutral excepted, by design).
stopifnot(all(semantic %in% ATTEND_PALETTE))

## --- the highlight is ONE group, and it resolves across id spaces ------------------
## attend_classes.R is not base-R sourceable (it loads tidyverse), so the config is read as
## text — the same trick test_plot_style.R rule [5] uses.
cls <- readLines(file.path("code", "attend_classes.R"), warn = FALSE)
hl  <- grep("^attend_highlight\\s*<-", cls)
stopifnot(length(hl) == 1L)
hl_block <- cls[hl:min(hl + 25L, length(cls))]
hl_code  <- hl_block[!grepl("^\\s*#", hl_block)]

## The retired `cohort` group must not come back as a highlight. It covered ~16 of ~40
## patients — a second competing fill rather than a mark — and encoded data availability,
## which the master already carries as the in_* columns and report 01 already draws.
stopifnot(!any(grepl("cohort\\s*=\\s*list", hl_code)))
stopifnot(any(grepl('polipo\\s*=\\s*list', hl_code)),
          any(grepl(HIGHLIGHT_BLUE, hl_code, fixed = TRUE)))

## No report may switch the overlay off. Marking one sample is only useful if the reader can
## find it in EVERY distribution it appears in; reports 05 and 06 each carried a wrapper that
## silently dropped it, which is indistinguishable on the page from never having had it.
for (f in list.files(file.path("analysis"), pattern = "\\.Rmd$", full.names = TRUE)) {
  txt <- readLines(f, warn = FALSE)
  bad <- grep("highlight_points\\s*<-\\s*function|highlight\\s*=\\s*NULL", txt)
  bad <- bad[!grepl("^\\s*#", txt[bad])]
  if (length(bad))
    stop(basename(f), ":", paste(bad, collapse = ","),
         ": this report overrides or disables highlight_points(). The polipo overlay is drawn ",
         "on every per-patient figure; there is no per-report opt-out.")
}

## 21S188 is a BARCODE and the master is keyed on pid, so a literal match found nothing and
## the overlay drew nothing without erroring. highlight_expand_ids() closes that gap.
stopifnot(is.function(register_highlight_xwalk), is.function(highlight_expand_ids))
register_highlight_xwalk(data.frame(barcode = c("21S188", "19S30"),
                                    pid     = c("P021", "P019"),
                                    stringsAsFactors = FALSE))
stopifnot("P021" %in% highlight_expand_ids("21S188"),   # barcode -> pid
          "21S188" %in% highlight_expand_ids("P021"),   # and back
          !"P019" %in% highlight_expand_ids("21S188"))  # never leaks across patients

hl_cfg <- list(groups = list(polipo = list(ids = "21S188", color = HIGHLIGHT_BLUE)))
got <- highlight_group_of(data.frame(pid = c("P021", "P019", "P007"),
                                     stringsAsFactors = FALSE), highlight = hl_cfg)
stopifnot(identical(got, c("polipo", NA_character_, NA_character_)))
stopifnot(grepl("1/1 matched", highlight_coverage(data.frame(pid = "P021",
                                                            stringsAsFactors = FALSE),
                                                 highlight = hl_cfg)))
## Empty registry -> literal matching only, exactly the pre-crosswalk behaviour.
register_highlight_xwalk()
stopifnot(identical(highlight_expand_ids("21S188"), "21S188"))
cat("palette + highlight: verified\n")

## --- attend_highlight_show: switching the MARK off must not drop the PATIENT ---------
## The flag is ONE global in attend_classes.R, not a per-report argument: the ban above is
## on a report opting out alone, which left the mark present in some distributions and
## absent in others. A global keeps every report saying the same thing either way.
hl_flag <- grep("^attend_highlight_show\\s*<-", cls, value = TRUE)
if (length(hl_flag) != 1)
  stop("attend_classes.R must define attend_highlight_show exactly once (found ",
       length(hl_flag), ") — the highlight switch is global, never per-report.")
stopifnot(any(grepl("^attend_highlight_show\\s*<-\\s*(TRUE|FALSE)\\s*$", hl_flag)))

## A CALIBRATION check, not a unit check. highlight_points() splits the cohort into the
## highlighted rows and the rest, and `base` is the is.na(.hl) half — so gating the mark
## anywhere after that split silently DELETES the highlighted patient from the figure. At
## every composition call site attend_box() is passed points = FALSE and highlight_points()
## IS the points layer, so the loss is a missing patient, drawn without an error. That is
## the same failure the id-expansion block above exists for: an overlay that fails quietly.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  register_highlight_xwalk(data.frame(barcode = "21S188", pid = "P021",
                                      stringsAsFactors = FALSE))
  d_hl <- data.frame(pid = c("P021", "P019", "P007", "P003"), g = c("a", "a", "b", "b"),
                     y = c(1, 2, 3, 4), stringsAsFactors = FALSE)
  plotted <- function(show) {
    attend_highlight_show <<- show
    b <- ggplot2::ggplot_build(ggplot2::ggplot(d_hl, ggplot2::aes(g, y)) +
                                 highlight_points(d_hl, "g", "y", highlight = hl_cfg))
    sum(vapply(b$data, nrow, integer(1)))
  }
  on  <- plotted(TRUE)
  off <- plotted(FALSE)
  if (on != nrow(d_hl) || off != nrow(d_hl))
    stop("attend_highlight_show drops patients: ", on, " plotted with the mark on, ",
         off, " with it off, of ", nrow(d_hl), ". The flag suppresses the MARK only.")
  rm(attend_highlight_show)
  register_highlight_xwalk()
  cat("highlight switch: mark toggles, cohort does not (", on, "/", off, " of ",
      nrow(d_hl), ")\n", sep = "")
}

## --- attend_box(): the mark is chosen per group, from the data ---------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  d <- data.frame(g = rep(c("big", "small"), times = c(20, 4)),
                  v = c(seq_len(20), 1:4), stringsAsFactors = FALSE)
  cls <- vapply(attend_box(d, x = "g", y = "v"), function(l) class(l$geom)[1], character(1))
  stopifnot("GeomBoxplot" %in% cls,    # the n = 20 group keeps its quartiles
            "GeomCrossbar" %in% cls,   # the n = 4 group gets a median bar instead
            "GeomPoint" %in% cls,      # points are ALWAYS drawn
            "GeomText" %in% cls)       # and n is always stated

  # points = FALSE for the call sites that supply highlight_points() themselves.
  cls2 <- vapply(attend_box(d, x = "g", y = "v", points = FALSE),
                 function(l) class(l$geom)[1], character(1))
  stopifnot(!"GeomPoint" %in% cls2)

  # All-large data must not emit a crossbar; all-small must not emit a box.
  big <- data.frame(g = rep(c("a","b"), each = 15), v = rnorm(30))
  sml <- data.frame(g = rep(c("a","b"), each = 3),  v = rnorm(6))
  cb <- vapply(attend_box(big, "g", "v"), function(l) class(l$geom)[1], character(1))
  cs <- vapply(attend_box(sml, "g", "v"), function(l) class(l$geom)[1], character(1))
  stopifnot(!"GeomCrossbar" %in% cb, !"GeomBoxplot" %in% cs)
  cat("attend_box: per-group mark selection verified\n")
}

## --- attend_theme() returns a ggplot theme (only if ggplot2 is installed) ---
if (requireNamespace("ggplot2", quietly = TRUE)) {
  th <- attend_theme()
  stopifnot(inherits(th, "theme"), inherits(th, "gg"))
} else {
  cat("attend_theme: ggplot2 absent — theme construction not exercised locally\n")
}

## --- SCNA plotters still build without error after the theme/palette swap (Task 2) ---
## (Colour-correctness is asserted by the palette-constant test in Task 1 + the grep
## guard in Step 4; ggplot scale introspection is too version-brittle to assert here,
## so this is a construction smoke test only.)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  pk <- data.frame(peak_id = "p1", chrom = "chr1",
                   wide_start = 1e6, wide_end = 2e6, direction = "amp",
                   stringsAsFactors = FALSE)
  fq <- data.frame(peak_id = rep("p1", 4),
                   scna_group = c("MMRd-high","MMRp-high","MMRd-low","MMRp-low"),
                   freq = c(0.6, 0.2, 0.1, 0.05), stringsAsFactors = FALSE)
  g <- scna_delta_plot(fq, pk, title = "t")
  stopifnot(inherits(g, "ggplot"))
  cat("scna_delta_plot: builds after palette/theme swap\n")
}

## --- attend_fig_save(): the 7pt swap must not eat the call site's theme tweaks ---
## theme_pub() is a COMPLETE ggplot theme, so `p + theme_pub()` REPLACES the accumulated
## theme. Without the carry-over in attend_fig_save(), a report figure exported for the
## manuscript silently regained the legend it had switched off and lost its rotated axis
## labels. Caught by rendering the export and looking at it; pinned here so it stays fixed.
if (requireNamespace("ggplot2", quietly = TRUE) &&
    file.exists(file.path("code", "plotstyle.R"))) {
  d <- data.frame(g = rep(c("a", "b"), each = 12), v = seq_len(24), stringsAsFactors = FALSE)
  p <- ggplot2::ggplot(d, ggplot2::aes(g, v, fill = g)) +
    ggplot2::geom_boxplot() + attend_theme() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
  stem <- file.path(tempdir(), "attend_fig_save_test")
  w <- attend_fig_save(p, stem, width = "single", data = d, formats = "png")
  stopifnot(file.exists(paste0(stem, ".png")),
            file.exists(paste0(stem, "_source_data.csv")),
            length(w) == 2L)

  # Rebuild what attend_fig_save() would have produced and inspect its theme directly.
  local({
    source(file.path("code", "plotstyle.R"), local = TRUE)
    th <- (p + theme_pub() +
             ggplot2::theme(legend.position = p$theme$legend.position,
                            axis.text.x = ggplot2::element_text(
                              angle = p$theme$axis.text.x$angle,
                              hjust = p$theme$axis.text.x$hjust)))$theme
    stopifnot(identical(th$legend.position, "none"),   # survived
              identical(th$axis.text.x$angle, 20))     # survived
    # ...and the point of exporting at all: the base size really did drop to print size.
    stopifnot(th$text$size == BASE_PT)
  })
  unlink(c(paste0(stem, ".png"), paste0(stem, "_source_data.csv")))
  cat("attend_fig_save: theme carry-over + 7pt swap verified\n")
}

## --- attend_aneu_cols carries the MMRd-subclass spellings (Phase 1.5) ---
stopifnot(attend_aneu_cols[["MMRd aneuploidy-low"]]  == attend_aneu_low,
          attend_aneu_cols[["MMRd aneuploidy-high"]] == attend_aneu_high)

cat("test_figure_system: ALL PASS\n")
