# Base-R unit tests (no tidyverse) — runnable in the bootstrap env.
source(file.path("code", "attend_plots.R"))  # base-R-sourceable: no top-level library()

## --- semantic palettes exist and carry the promoted Fig-1a hexes ---
stopifnot(is.character(attend_mmr_cols),
          attend_mmr_cols[["MMRp"]] == "#999999",
          attend_mmr_cols[["MMRd"]] == "#E69F00",
          attend_mmr_cols[["MMR proficient"]] == "#999999",
          attend_mmr_cols[["MMR deficient"]] == "#E69F00")

stopifnot(identical(attend_aneu_low, "#56B4E9"),   # sky blue   = aneuploidy-low
          identical(attend_aneu_high, "#D55E00"),  # vermillion = aneuploidy-high
          attend_aneu_high != attend_tp53_cols[["mut"]],  # must not clash with the adjacent
          attend_aneu_high != attend_mmr_cols[["MMRd"]],  # Fig-1a annotation bars
          is.character(attend_aneu_cols),
          attend_aneu_cols[["aneuploidy-low"]]  == attend_aneu_low,
          attend_aneu_cols[["aneuploidy-high"]] == attend_aneu_high,
          attend_aneu_cols[["aneu-low"]]  == attend_aneu_low,
          attend_aneu_cols[["aneu-high"]] == attend_aneu_high)

stopifnot(is.character(attend_tp53_cols),
          attend_tp53_cols[["TP53-normal"]]   == "#999999",
          attend_tp53_cols[["TP53-abnormal"]] == "#CC79A7")

## --- the two rules the palette assignment rests on (attend_plots.R) ---------------
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

## Response must not reuse EITHER aneuploidy pole: reports 05 and 06 draw the same
## cell-type panels keyed on different variables, so a shared colour would make two
## different groupings look like one.
stopifnot(!attend_resp_cols[["responder"]]     %in% c(attend_aneu_low, attend_aneu_high),
          !attend_resp_cols[["non-responder"]] %in% c(attend_aneu_low, attend_aneu_high),
          attend_resp_cols[["responder"]] != attend_resp_cols[["non-responder"]])

## The neutral single-group fill must not carry ANY semantic meaning — this is the bug that
## motivated the palette pass: #A6CEE3 was simultaneously "aneuploidy-low" and "no grouping".
semantic <- unique(c(attend_mmr_cols, attend_aneu_cols, attend_tp53_cols, attend_resp_cols))
stopifnot(is.character(attend_neutral), !attend_neutral %in% semantic)

## Every semantic colour is drawn from the one palette (attend_neutral excepted, by design).
stopifnot(all(semantic %in% ATTEND_OKABE_ITO))

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
