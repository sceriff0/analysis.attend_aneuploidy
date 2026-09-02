# ============================================================================
# code/attend_plots.R — TCGA Figure 1a/1b reproduction (Kandoth et al. 2013)
# Base-R-sourceable (NO top-level library()); heavy deps are requireNamespace-
# guarded inside the functions. Reports source() this AFTER attend_classes.R.
# ============================================================================

# Local null-coalesce so this file works under R 4.3.2 (base %||% arrived in 4.4)
# and standalone under base R. Harmless if tidyverse later redefines it.
`%||%` <- function(a, b) if (is.null(a)) b else a

# The canonical genome order for chromosome factors used by every Fig-1a heatmap.
.fig1a_chrom_levels <- paste0("chr", c(as.character(1:22), "X", "Y"))

# Map CNV feature names to (chrom, coord) so a heatmap can put chromosomal
# location on the vertical axis (TCGA Fig. 1a). Handles genome-bin colnames
# ("chr1:0"), arm names ("1p"/"17q"), and GISTIC peak genes (needs `cytoband`,
# a named gene->band vector like c(MYC="8q24.21")). Returns a data.frame
# feature/chrom/coord ordered by genome position. Base R (no tidyverse) so it is
# unit-testable without the full stack.
feature_positions <- function(features, cytoband = NULL) {
  features <- as.character(features)
  chrom <- rep(NA_character_, length(features)); coord <- rep(NA_real_, length(features))

  is_bin <- grepl("^chr[0-9XY]+:", features, ignore.case = TRUE)
  is_arm <- grepl("^[0-9]{1,2}[pq]$|^[XY][pq]$", features)

  # bins: "chr1:1000000"
  if (any(is_bin)) {
    sp <- strsplit(features[is_bin], ":", fixed = TRUE)
    chrom[is_bin] <- vapply(sp, function(x) x[[1]], character(1))
    coord[is_bin] <- as.numeric(vapply(sp, function(x) x[[2]], character(1)))
  }
  # arms: "17q" -> chr17, coord p=1/q=2
  if (any(is_arm)) {
    a <- features[is_arm]
    chrom[is_arm] <- paste0("chr", sub("[pq]$", "", a))
    coord[is_arm] <- ifelse(grepl("p$", a), 1, 2)
  }
  # remaining = gene symbols resolved via cytoband ("8q24.21" -> chr8, band ordinal)
  rest <- is.na(chrom)
  if (any(rest)) {
    band <- if (!is.null(cytoband)) cytoband[features[rest]] else NA_character_
    m <- regmatches(band, regexec("^([0-9XY]{1,2})([pq])([0-9.]*)", band))
    chrom[rest] <- vapply(m, function(x) if (length(x) >= 2) paste0("chr", x[[2]]) else NA_character_, character(1))
    coord[rest] <- vapply(m, function(x) {
      if (length(x) < 4) return(NA_real_)
      arm_ord <- if (identical(x[[3]], "p")) 0 else 1000       # q sits after all p
      arm_ord + suppressWarnings(as.numeric(ifelse(nzchar(x[[4]]), x[[4]], "0")))
    }, numeric(1))
  }

  out <- data.frame(feature = features, chrom = chrom, coord = coord, stringsAsFactors = FALSE)
  out <- out[!is.na(out$chrom), , drop = FALSE]
  out$chrom <- factor(out$chrom, levels = .fig1a_chrom_levels)
  out <- out[!is.na(out$chrom), , drop = FALSE]
  out[order(out$chrom, out$coord), , drop = FALSE]
}

# ============================================================================
# THE COLOUR SYSTEM — one palette, one place.
# ============================================================================
# Everything categorical in this pipeline resolves to ONE palette: Paul Tol's
# "vibrant" qualitative scheme (Tol, "Colour Schemes", SRON/EPS technical note
# v3.2, Fig. 3). Six hues, distinguishable under all three common colour-vision
# deficiencies and separable in greyscale by lightness. Before this block existed
# the pipeline mixed six palettes (Okabe-Ito, ColorBrewer Paired, ColorBrewer
# Set1/RdBu, matplotlib tab10, seaborn deep, and survminer's named R colours),
# which is why the same light blue meant "aneuploidy-low" in report 05 and "no
# grouping at all" in reports 01/09.
#
# WHY TOL VIBRANT AND NOT OKABE-ITO. The two poles of aneuploidy have to read as
# LIGHT BLUE (quiet) and RED (altered) — that is the direction the SCNA heatmap
# body beneath the annotation bar already reads, and the pair the project uses in
# every figure legend and in the manuscript text. Okabe-Ito has no true red: its
# warm end is vermillion #D55E00, which reads as a dark orange next to MMRd's
# #E69F00 and forced the aneuploidy scale to be described as something it was not.
# Tol vibrant carries a genuine red (#CC3311) AND a genuine light blue (#33BBEE)
# while remaining CVD-safe as a set, so the honest colour and the accessibility
# constraint stop competing.
#
# Two deliberate deviations from Tol's published vector:
#   * Tol's grey is #BBBBBB. That is visually indistinguishable from
#     `attend_neutral` (#BFBFBF, the "no grouping at all" fill), which would
#     recreate the exact bug this palette pass fixed: a non-semantic box reading as
#     a semantic reference state. The reference grey therefore stays the darker
#     mid-grey #999999.
#   * `black` is appended as an 8th slot. It is NEVER a semantic colour — it exists
#     only so `attend_pal()` and the k = 2..8 cluster bar keep eight distinct
#     values without recycling.
#
# Two rules make the semantic assignments below safe:
#   [A] GREY IS ALWAYS THE REFERENCE STATE. MMRp, TP53-normal, wild-type — the
#       thing a reader should not look at — is grey in every figure. Only the
#       altered state carries a hue.
#   [B] CONCEPTS THAT CAN SHARE A FIGURE MUST NOT SHARE A HUE. The Fig-1a
#       annotation stack puts MMR, aneuploidy, TP53 and cluster side by side, so
#       those four are mutually distinct. TP53 and response reuse #EE3377 because
#       they provably never co-occur: report 06 is the only report that colours by
#       response and it contains no TP53, and no report that plots TP53 (00/05/07/09)
#       colours by response. Adding TP53 to report 06 means giving one of them a
#       new hue — see test_plot_style.R, which pins this.
#   [C] `blue` (#0077BB) IS RESERVED FOR THE HIGHLIGHT. It is the one hue no
#       semantic palette may claim: highlight points are drawn ON TOP of boxes
#       filled with the semantic colours, so a highlight sharing a hex is invisible
#       on exactly the box it exists to mark. Keeping one hue unclaimed is what
#       makes that guarantee cheap. test_figure_system.R pins it.
# `pale_blue`/`salmon` are not Tol's: they are ColorBrewer #A6CEE3 / #E31A1C as they render
# at ~50% alpha, carried in the palette rather than written at the call site so the rule
# "colour comes from the palette, never from the call site" still holds. `cyan` and `red`
# stay in the vector — attend_pal() and the cluster bar draw from it positionally, and
# removing two slots would resequence every cluster colour in reports 09 and 10.
ATTEND_PALETTE <- c(blue      = "#0077BB",   # RESERVED — the polipo highlight, rule [C]
                    cyan      = "#33BBEE",   # (was aneuploidy-low; kept for attend_pal())
                    teal      = "#009988",   # responder
                    orange    = "#EE7733",   # MMRd
                    red       = "#CC3311",   # (was aneuploidy-high; kept for attend_pal())
                    magenta   = "#EE3377",   # TP53-abnormal, and non-responder (rule [B])
                    grey      = "#999999",   # every reference state, rule [A]
                    black     = "#000000",   # non-semantic filler, cluster bars only
                    pale_blue = "#D6E5F0",   # aneuploidy-low
                    salmon    = "#E3918F")   # aneuploidy-high

# TCGA Fig. 1a colour for the copy-number clusters, matching the paper's Cluster bar.
# Covers up to 8 clusters so the k = 2..8 sweep (fig1a_heatmap_ksweep()) always has a
# distinct colour per cluster; the standard pipeline fixes k at 4 (Fig. 2b). k > 8
# would need more colours. Derived from ATTEND_PALETTE rather than the tab10 set it
# used to hardcode, so a cluster bar and a ggplot legend in the same report agree.
.fig1a_cluster_cols <- stats::setNames(unname(ATTEND_PALETTE), as.character(1:8))

# Project-wide qualitative palette — the brewer-free replacement for
# scale_*_brewer(palette = "Set1") / RColorBrewer::brewer.pal(). Reports use
#   scale_colour_manual(values = attend_pal(n), aesthetics = c("colour", "fill"))
# instead, so no plot depends on RColorBrewer being loadable.
# `n` > 8 recycles with a warning — a categorical plot needing more than 8 colours
# should use a different encoding anyway. Pass `names` for a named vector.
attend_pal <- function(n = 8, names = NULL) {
  if (!is.null(names)) n <- length(names)
  if (n > length(ATTEND_PALETTE))
    warning("attend_pal(): ", n, " colours requested but only ", length(ATTEND_PALETTE),
            " distinct ones exist — recycling.", call. = FALSE)
  out <- unname(rep(ATTEND_PALETTE, length.out = max(1, n)))
  if (!is.null(names)) stats::setNames(out, names) else out
}

# --- Semantic colour palettes (single source of truth) --------------------
# Keyed by every label spelling in use across reports so one manual scale covers the
# Fig-1a annotation bars and the ggplots alike.

# Reference state, rule [A]. Also the fill for a boxplot with no grouping at all:
# a single-group box must be visibly NON-semantic, or a reader carries a meaning
# into it from the last figure. This is deliberately not any hue above.
attend_neutral <- "#BFBFBF"

attend_mmr_cols <- c(MMRp = unname(ATTEND_PALETTE["grey"]),
                     MMRd = unname(ATTEND_PALETTE["orange"]),
                     `MMR proficient` = unname(ATTEND_PALETTE["grey"]),
                     `MMR deficient`  = unname(ATTEND_PALETTE["orange"]))

# Aneuploidy is a two-pole scale, not a reference/altered pair, so both poles carry a
# hue. LIGHT BLUE (quiet) -> RED (altered) reads the same direction as the SCNA heatmap
# body beneath the annotation bar, so the bar and the matrix do not tell opposite colour
# stories. This is the pair the project describes in prose and in the manuscript, and
# Tol vibrant is the reason it can now be literal: #33BBEE is a light blue and #CC3311
# is a red, where the earlier Okabe-Ito stand-ins (#56B4E9 sky blue / #D55E00
# vermillion) forced the warm pole to be called red while rendering as dark orange next
# to MMRd's own orange.
# The pre-f77f7e8 pair, restored on request so the composition figures match the reference
# figure the manuscript is built around. Sampled from it: these are ColorBrewer #A6CEE3 /
# #E31A1C as they RENDER at ~50% alpha. Baked in at full opacity rather than set as an alpha,
# so the points and the polipo highlight drawn on top stay crisp instead of compositing.
#
# The documented reason #A6CEE3 was retired -- it was simultaneously "aneuploidy-low" and
# "no grouping", so unrelated boxplots read as aneuploidy-low -- no longer applies: the
# no-grouping fill is attend_neutral #BFBFBF and has been since the neutral was split out.
# The light-blue (quiet) -> red (altered) direction the repalette existed to protect is
# unchanged, and both stay separable in greyscale (luminance 228 vs 168) and against MMRd's
# orange. Rule [C] holds: the highlight blue #0077BB is solid and larger, so it reads on top
# of the pale box rather than into it.
attend_aneu_low  <- unname(ATTEND_PALETTE["pale_blue"])
attend_aneu_high <- unname(ATTEND_PALETTE["salmon"])
attend_aneu_cols <- c(`aneuploidy-low` = attend_aneu_low, `aneuploidy-high` = attend_aneu_high,
                      `aneu-low` = attend_aneu_low, `aneu-high` = attend_aneu_high,
                      `MMRd aneuploidy-low` = attend_aneu_low,
                      `MMRd aneuploidy-high` = attend_aneu_high)

# TP53: reference/altered pair, so grey + one hue (rule [A]). Magenta is the only warm
# hue left once aneuploidy has taken red and MMR has taken orange, and it is the one
# that stays separable from both in the adjacent Fig-1a bars — #EE3377 differs from
# #CC3311 in hue AND from #EE7733 in hue, so the three-bar stack does not read as one
# warm smear under deuteranopia.
attend_tp53_cols <- c(`TP53-normal` = unname(ATTEND_PALETTE["grey"]),
                      `TP53-abnormal` = unname(ATTEND_PALETTE["magenta"]),
                      wt = unname(ATTEND_PALETTE["grey"]),
                      mut = unname(ATTEND_PALETTE["magenta"]))

# Treatment response (add_response_class). Response is plotted on the SAME cell-type
# panels as aneuploidy in reports 05/06, so it must not reuse the aneuploidy pair —
# two different groupings would look like one grouping at a glance.
attend_resp_cols <- c(`responder` = unname(ATTEND_PALETTE["teal"]),
                      `non-responder` = unname(ATTEND_PALETTE["magenta"]),
                      `MMRd responder` = unname(ATTEND_PALETTE["teal"]),
                      `MMRd non-responder` = unname(ATTEND_PALETTE["magenta"]))

# Cohort of origin (report 11's ATTEND-vs-TCGA contrast). This is rule [A] applied to a
# non-biological variable: TCGA UCEC is the EXTERNAL REFERENCE cohort — the thing the ATTEND
# result is measured against — so it takes the reference grey, and only ATTEND carries a hue.
# Teal, because red and light blue are already the aneuploidy poles that report 11 plots
# alongside this, and orange is MMR (report 11 facets on MMR group). Teal is otherwise
# claimed only by `responder`, and no report colours by cohort and by response at once.
attend_cohort_cols <- c(`TCGA UCEC (primary)`  = unname(ATTEND_PALETTE["grey"]),
                        `ATTEND (metastatic)`  = unname(ATTEND_PALETTE["teal"]),
                        TCGA   = unname(ATTEND_PALETTE["grey"]),
                        ATTEND = unname(ATTEND_PALETTE["teal"]))

# ============================================================================
# THE THEME — one base_size, one grid policy.
# ============================================================================
# Reports render to a workflowr HTML site, so this is sized for a screen (11pt), not
# for a journal column. The manuscript path is attend_fig_save(), which restyles to
# plotstyle.R's 7pt theme_pub() on the way out — see below.
#
# Grid policy: horizontal major only. Every figure in this pipeline reads a VALUE off
# the y axis (a fraction, a score, mut/Mb); none reads one off a categorical x, so an
# x gridline is ink that guides nothing. Minor gridlines are off for the same reason.
attend_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = ggplot2::element_text(size = base_size - 1, colour = "grey30"),
      plot.caption     = ggplot2::element_text(size = base_size - 2, colour = "grey30", hjust = 0),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size - 1),
      strip.background = ggplot2::element_blank(),
      axis.title       = ggplot2::element_text(size = base_size - 1),
      axis.text        = ggplot2::element_text(size = base_size - 2, colour = "grey20"),
      legend.title     = ggplot2::element_text(size = base_size - 1),
      legend.text      = ggplot2::element_text(size = base_size - 2),
      legend.key.size  = ggplot2::unit(4, "mm"),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.3, colour = "grey88"))
}

# ============================================================================
# attend_box() — a distribution layer that tells the truth about small n.
# ============================================================================
# ATTEND is a ~40-patient cohort; split by MMR, by aneuploidy class, and by responder
# status, a group of four is routine. A boxplot of four points draws quartiles that are
# an artefact of which two points landed in the middle, and `outlier.shape = NA` — which
# every call site in this pipeline passes, to avoid double-plotting points — DELETES the
# outliers from a group too small to have any.
#
# So the mark is chosen per group from the data, at knit time:
#   n >= min_n  ->  box + jittered points   (quartiles are supported)
#   n <  min_n  ->  jittered points + a median crossbar, no box
# Both branches always draw the points, so nothing is ever summarised away invisibly.
#
# `by` carries the FACET variables. Group size differs per panel, so a cell type with
# 20 patients in one facet and 4 in another must get the box in the first and the
# crossbar in the second; keying the split on x alone would give both the same verdict.
#
# Returns a list of ggplot layers, so it composes: ggplot(d, aes(...)) + attend_box(...).
# Base R only (aggregate, not dplyr) — attend_plots.R is sourced in the bootstrap env.
# 5, not 10. ATTEND's aneuploidy-high arm is ~5 patients, and at 10 that arm drew as a bare
# median line beside a full box — which reads as "no data here" rather than "few patients
# here", and made every composition panel look collapsed. Quartiles from 5 points are noisy
# and that is a real cost, but the n is printed under every group, so the reader can discount
# a box the same way they would discount the crossbar. Below 5 the crossbar still applies.
attend_min_box_n <- 5L

# `linewidth` and `label_size` are arguments rather than constants because a theme cannot
# reach either: the box border is a geom linewidth in mm and the n= label is a geom_text
# size in mm, so attend_theme(base_size =) scales the axes around them and leaves both at
# their site defaults. The composition figures pass ihc_fig_style()'s values.
attend_box <- function(data, x, y, by = NULL, min_n = attend_min_box_n,
                       jitter_width = 0.18, point_size = 0.8, point_alpha = 0.55,
                       box_width = 0.6, label_n = TRUE, point_colour = "black",
                       fill = attend_neutral, points = TRUE, expand_y = label_n,
                       linewidth = 0.5, label_size = 2.4) {
  keys <- c(x, by)
  d <- data[stats::complete.cases(data[, c(keys, y), drop = FALSE]), , drop = FALSE]
  if (nrow(d) == 0) return(list())

  # n per (x, facet...) cell, then the two row-sets the layers draw from.
  cnt <- stats::aggregate(d[[y]], by = as.list(d[, keys, drop = FALSE]),
                          FUN = function(v) sum(is.finite(v)))
  names(cnt)[ncol(cnt)] <- ".n"
  cell <- function(df) do.call(paste, c(lapply(keys, function(k) as.character(df[[k]])), sep = "\r"))
  dense_cells  <- cell(cnt)[cnt$.n >= min_n]
  d_dense  <- d[cell(d) %in% dense_cells, , drop = FALSE]
  d_sparse <- d[!cell(d) %in% dense_cells, , drop = FALSE]

  layers <- list()
  if (nrow(d_dense) > 0) {
    # fill = NULL means "inherit the parent aes(fill = ...)" — for the rare panel whose
    # fill encodes something the x axis does NOT already carry. Every other call gets the
    # neutral grey, because a box coloured by its own x position is redundant ink.
    box_args <- list(data = d_dense, outlier.shape = NA, width = box_width,
                     linewidth = linewidth)
    if (!is.null(fill)) box_args$fill <- fill
    layers <- c(layers, list(do.call(ggplot2::geom_boxplot, box_args)))
  }
  if (nrow(d_sparse) > 0)
    # A median crossbar, drawn only where the box was withheld. Same width as the box
    # so the two branches line up when they appear in adjacent facets.
  {
    cb_args <- list(data = d_sparse, fun = stats::median, geom = "crossbar",
                    width = box_width, linewidth = linewidth)
    if (!is.null(fill)) cb_args$fill <- fill
    layers <- c(layers, list(do.call(ggplot2::stat_summary, cb_args)))
  }

  # points = FALSE for the call sites that supply their own overlay — highlight_points()
  # partitions the cohort into highlight groups and a base layer, and drawing this jitter
  # underneath it would plot every patient twice.
  if (points)
    layers <- c(layers, list(ggplot2::geom_jitter(
      width = jitter_width, height = 0, size = point_size, alpha = point_alpha,
      colour = point_colour, show.legend = FALSE)))

  if (label_n) {
    # n at the panel floor: y = -Inf keeps it anchored whatever the free_y range is.
    cnt$.lab <- paste0("n=", cnt$.n)
    layers <- c(layers, list(ggplot2::geom_text(
      data = cnt, inherit.aes = FALSE,
      mapping = ggplot2::aes(x = .data[[x]], y = -Inf, label = .data$.lab),
      vjust = -0.6, size = label_size, colour = "grey35")))
    # ggplot's default 5% bottom expansion is not enough room for a label at y = -Inf: a
    # group whose points sit on the panel floor collides with its own n. Widen the bottom
    # only. Pass expand_y = FALSE at a call site that sets its own y scale (e.g.
    # scale_y_log10) — the later scale wins anyway, but silently, with a replacement message.
    if (expand_y)
      layers <- c(layers, list(ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0.12, 0.05)))))
  }
  layers
}

# ============================================================================
# shared_frac_y() / frac_display_y() — one axis for every facet in a figure.
# ============================================================================
# The IHC cell-type panels faceted on `scales = "free_y"`, which rescales every panel to its
# own data: a cell type at a fraction of a percent of the tissue and one at half the tissue
# drew the same box in the same place, so the panel RANGE carried the abundance and the
# reader could not see it. One shared axis over every cell type in the figure fixes that, and
# is what the reports ask for: Macrophages must be readable against T cytotoxic on the same
# scale, not against a rescaled axis of its own.
#
# An intermediate design split the facets into abundance tiers and gave each tier its own
# figure and its own shared axis. It is gone — a reader comparing two cell types should not
# have to check which figure each landed in, and the tier boundaries were a property of the
# cohort rather than of the biology. The cost is real and is accepted deliberately: on the
# tumour-cell denominator these cell types span orders of magnitude, so the rarest panels sit
# low against a scale set by the most abundant. frac_display_y() is what keeps that from
# being hopeless — the axis stops at the largest whisker, not the largest point.
#
# The p-value label MUST be placed explicitly, by frac_display_y(). ggpubr derives its
# default `label.y` from the PLOT-WIDE y range, not from the panel's own data, so on a shared
# axis every facet gets the same default, pinned to the global maximum. Measured on a four-
# facet repro with one patient at 10x the rest: all four labels land at y = 0.60 while their
# boxes sit under 0.09 — adrift at the ceiling, which on the page reads as no p-value at all.
# An earlier fix deleted an explicit label_y on the grounds that ggpubr "places each label
# just above the data it belongs to". That is true under `scales = "free_y"`, which these
# facets used until the shared axis landed, and false after. The label height and the axis
# ceiling are therefore ONE decision, and frac_display_y() returns both from one call.
#
# frac_display_y() also removes what set the labels adrift: a single extreme patient ranging
# the whole figure. `outlier.shape = NA` hides the outlier glyph but not its pull on the
# scale, and highlight_points() redraws every patient regardless. The ceiling is the largest
# upper whisker over all (x, facet) cells, so no box and no whisker is ever cut, and it is
# applied with coord_cartesian() — a ZOOM. scale_y_continuous(limits=) would drop the
# off-scale rows BEFORE stat_boxplot computes and silently move the medians it is supposed to
# be showing. Points above the ceiling are counted, never quietly dropped: `$note` carries
# the count into the subtitle.
#
# attend_box(label_n = TRUE) appends its own scale_y_continuous() for the n= labels at
# y = -Inf, and a second one REPLACES it, silently — so call sites pass `expand_y = FALSE`
# and take the bottom expansion from here.
#
# Base R only — attend_plots.R is sourced in the bootstrap env, so ggplot2 may only be
# touched at call time.
# The display range for one shared-axis figure, and the label height that has to sit inside
# it. `keys` are the columns whose combination is one box (the x variable and the facet).
# Returns $ylim, $label_y, $n_above and a $note for the subtitle; every field is NULL/0/""
# when there is nothing finite to scale, so a call site can pass the result through blind.
frac_display_y <- function(df, value_col, keys, headroom = 0.18, whisker_mult = 1.5) {
  none <- list(ylim = NULL, label_y = NULL, n_above = 0L, note = "")
  v  <- suppressWarnings(as.numeric(df[[value_col]]))
  ok <- is.finite(v)
  if (!any(ok)) return(none)

  # One upper whisker per (x, facet) cell. Same \r-joined key idiom as attend_box(), which
  # has to survive factor levels that contain the separator.
  cells <- do.call(paste, c(lapply(keys, function(k) as.character(df[[k]])[ok]), sep = "\r"))
  uw <- tapply(v[ok], cells, function(x) {
    q <- stats::quantile(x, c(0.25, 0.75), names = FALSE, na.rm = TRUE)
    q[2] + whisker_mult * (q[2] - q[1])
  })
  uw  <- uw[is.finite(uw)]
  top <- if (length(uw)) max(uw) else max(v[ok])
  lo  <- min(0, min(v[ok]))
  # Nothing above the floor (every value at or below 0) leaves no range to zoom into.
  # A flat set at a non-zero level is NOT degenerate: [0, v*(1+headroom)] is the right
  # window for it, and the boxes land where the reader expects.
  if (!is.finite(top) || top <= lo) return(none)

  span <- top - lo
  ceil <- lo + span * (1 + headroom)
  n_above <- sum(v[ok] > ceil)
  list(ylim    = c(lo, ceil),
       label_y = lo + span * (1 + headroom * 0.45),
       n_above = n_above,
       # Terse on purpose: this rides in the subtitle, and the long form cost two extra
       # wrapped lines. The guarantee it abbreviates -- nothing is dropped, the zoom is
       # coord_cartesian() -- is stated once in the derivation note above the chunk.
       note    = if (n_above > 0)
                   sprintf("; %d point%s off scale (all still in the boxes and the test)",
                           n_above, if (n_above == 1L) "" else "s")
                 else "")
}

# Pass the frac_display_y() result to apply its zoom; called with no argument this is the
# bare expansion it has always been.
shared_frac_y <- function(y = NULL, bottom = 0.12, top = 0.10) {
  out <- list(ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(bottom, top))))
  if (!is.null(y) && !is.null(y$ylim))
    out <- c(out, list(ggplot2::coord_cartesian(ylim = y$ylim)))
  out
}

# ============================================================================
# ihc_fig_style() — one scale for every size in a composition figure.
# ============================================================================
# The composition panels are read at a glance, so they carry larger type and heavier box
# borders than the site default. Four sizes have to grow together, and only ONE of them is
# owned by the theme:
#
#   base_size    -> attend_theme(), which sizes title, axes, strip text
#   label_size   -> geom_text size (mm) for the n= labels, inside attend_box()
#   p_size       -> stat_compare_means() size (mm), at the call site
#   linewidth    -> geom_boxplot/crossbar border (mm), inside attend_box()
#
# Scaling the theme alone leaves the n=, the p-value and the border at their site defaults,
# so the figure grows unevenly and the borders look thinner as the type gets bigger. Hence
# one bundle, derived from one number, rather than a multiplication at four call sites —
# the same reason frac_display_y() returns the ceiling and the label height together.
#
# Base R only, and it holds no ggplot objects: attend_plots.R is sourced in the bootstrap
# env, so ggplot2 may only be touched at call time.
# Defined HERE, not in attend_classes.R, and this is the k_cnv rule: attend_plots.R is
# sourced in the bootstrap env, so a default argument reaching into attend_classes.R dies
# with "object 'attend_ihc_text_scale' not found" for any caller that sources the plots
# layer alone -- every test in code/tests/ among them. Visibility follows the source()
# graph, not the directory listing. It is a figure constant, so it belongs beside
# attend_min_box_n.
attend_ihc_text_scale <- 2.5

# Wrap a figure title/subtitle to the canvas it will be drawn on. ggplot does NOT wrap and
# does NOT warn: an over-long title is clipped at the device edge, so at 2.5x type the first
# render lost the right-hand half of every title silently. `width_in` is the chunk's
# fig.width; 0.5 * size_pt is the mean glyph advance for this face, close enough that the
# result never overruns and rarely wraps a line early.
wrap_fig_text <- function(txt, width_in, size_pt, gutter = 0) {
  if (!length(txt) || !nzchar(txt)) return(txt)
  # 0.62, not 0.5: measured against the rendered title, a half-em advance over-estimates how
  # much fits and the line still ran off the canvas. `gutter` is the width the y-axis title
  # and tick labels take before the panel starts -- the plot title is left-aligned to the
  # PANEL, not the device, so a wrap computed on the full figure width overruns by exactly
  # that much. Callers pass it for titles and leave it 0 for an axis label.
  ch <- max(12L, as.integer(max(width_in - gutter, 1) * 72 / (0.62 * size_pt)))
  paste(strwrap(txt, width = ch), collapse = "\n")
}

# Tick labels for a two-group x axis carrying large type.
#
# At 2.5x, "aneuploidy-high" angled at 20 degrees runs off the LEFT edge of the device and
# ggplot clips it without warning -- the reported symptom was the leading "aneu" missing from
# the first label of every figure. Angling is also expensive vertically at this size, and
# that height comes straight out of the panels.
#
# So the labels are shortened and drawn HORIZONTAL. The prefix every level shares carries no
# information the title does not already give ("... by aneuploidy"), so
# aneuploidy-high/aneuploidy-low become high/low. Levels with nothing in common are left
# alone and merely broken over two lines: responder/non-responder keep both words, because
# there the words ARE the contrast.
short_class_labels <- function(x) {
  v <- as.character(x)
  if (length(v) > 1L && !anyNA(v) && all(nzchar(v))) {
    n <- min(nchar(v)); k <- 0L
    while (k < n && length(unique(substr(v, 1L, k + 1L))) == 1L) k <- k + 1L
    pre <- substr(v[1], 1L, k)
    # Trim the shared prefix back to its last separator, so a partial word is never cut.
    cut <- suppressWarnings(max(c(0L, gregexpr("[- ]", pre)[[1]])))
    if (cut > 0L && all(nchar(v) > cut)) v <- substring(v, cut + 1L)
  }
  ifelse(nchar(v) > 8L, sub("([- ])", "\\1\n", v), v)
}

ihc_fig_style <- function(scale = attend_ihc_text_scale, title_scale = 1) {
  for (nm in c("scale", "title_scale")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0)
      stop("ihc_fig_style(): `", nm, "` must be one positive finite number.")
  }
  list(scale      = scale,
       base_size  = 11  * scale,   # attend_theme()'s default
       label_size = 2.4 * scale,   # attend_box()'s n= label
       p_size     = 3   * scale,   # stat_compare_means()
       # Title and subtitle do NOT take `scale`. They are labelling, not data: at 2.5x they
       # cost three wrapped lines and a third of the canvas while telling the reader nothing
       # the strips and axis do not. They stay at the site's own size, so the builders must
       # override attend_theme(base_size=), which would otherwise scale them with everything
       # else. `title_scale` is the knob if they ever need to grow.
       title_size = 11 * title_scale + 1,
       sub_size   = 11 * title_scale - 1,
       # Width the y-axis title + tick labels consume before the panel begins. Grows with
       # the type, which is why it cannot be a constant at the call site.
       gutter_in  = 0.10 * (11 * scale - 1),
       # The border grows on a SQUARE ROOT of the scale, not linearly. At 2.5x a linear
       # 1.25mm border on a box a few mm tall reads as a filled slab and swallows the
       # median line; sqrt keeps it visibly heavier while the box still reads as a box.
       linewidth  = 0.5 * sqrt(scale))
}

# ============================================================================
# attend_fig_save() — the manuscript path.
# ============================================================================
# attend_theme() is sized for the HTML site. A figure headed for the paper needs
# plotstyle.R's 7pt theme_pub() at a real journal column width, plus the source-data
# CSV that makes the figure auditable. This restyles and hands off, so a report can
# render a figure to the site AND export it, from one plot object:
#
#   p <- ggplot(...) + attend_theme()
#   print(p)                                        # -> the site, 11pt
#   attend_fig_save(p, "figures/fig2_aneu", width = "double", data = plotted)
#
# plotstyle.R is only sourced when this is called, so the base-R bootstrap env that
# runs the unit tests never needs ggplot2 loaded at source() time.
attend_fig_save <- function(plot, path, width = "single", data = NULL, restyle = TRUE, ...) {
  ps <- Filter(file.exists, c(
    if (requireNamespace("here", quietly = TRUE)) here::here("code", "plotstyle.R"),
    file.path("code", "plotstyle.R"), "plotstyle.R"))
  if (length(ps) == 0)
    stop("attend_fig_save(): code/plotstyle.R not found — the manuscript style is missing.")
  source(ps[[1]], local = TRUE)
  if (restyle && inherits(plot, "ggplot")) {
    # theme_pub() is a COMPLETE theme, so `plot + theme_pub()` REPLACES the accumulated
    # theme rather than merging into it — which silently drops the call-site tweaks every
    # report adds on top of attend_theme(): legend.position = "none" and the rotated x-axis
    # labels that keep long class names legible. So they are read off the incoming plot and
    # re-applied after the swap. Only layout is carried over, never sizes: the whole point
    # of the export is to move from the site's 11pt to the journal's 7pt.
    keep <- list()
    th <- plot$theme
    if (!is.null(th$legend.position))  keep$legend.position  <- th$legend.position
    if (!is.null(th$legend.direction)) keep$legend.direction <- th$legend.direction
    for (ax in c("axis.text.x", "axis.text.y")) {
      e <- th[[ax]]
      if (inherits(e, "element_text") && (!is.null(e$angle) || !is.null(e$hjust)))
        keep[[ax]] <- ggplot2::element_text(angle = e$angle, hjust = e$hjust)
      if (inherits(e, "element_blank")) keep[[ax]] <- ggplot2::element_blank()
    }
    plot <- plot + theme_pub()
    if (length(keep)) plot <- plot + do.call(ggplot2::theme, keep)
  }
  pub_save(plot, path, width = width, data = data, ...)
}

# Threshold that splits a per-sample aneuploidy vector into high/low, used by BOTH the
# heatmap's binary aneuploidy bar and the aneuploidy boxplots' point shapes so the two agree.
# On the [0,1] master-score scale it uses the pipeline cutoff attend_thresholds$aneuploidy
# (matching aneuploidy_class everywhere else); on any other scale (e.g. ASCETS arm counts) it
# falls back to the cohort median. Returns NA when the vector has no finite values.
.aneu_split_threshold <- function(v) {
  v   <- suppressWarnings(as.numeric(v))
  rng <- range(v, na.rm = TRUE)
  if (!all(is.finite(rng))) return(NA_real_)
  on_unit <- rng[1] >= 0 && rng[2] <= 1
  if (on_unit && exists("attend_thresholds") && !is.null(attend_thresholds$aneuploidy))
    attend_thresholds$aneuploidy
  else stats::median(v, na.rm = TRUE)
}

# Readable, colour-blind-aware palettes for the ATTEND covariate bars drawn ON TOP of the
# Fig-1a heatmap. Keyed by annotation NAME; each entry maps the bar's factor levels (or a
# continuous ramp) to colours. This replaces ComplexHeatmap's auto-assigned defaults, which
# were low-contrast and drifted between reports. Colours are from the Okabe-Ito set so the
# categorical bars stay distinguishable in the common colour-vision deficiencies.

# The explicit label given to samples with no value, so "missing" appears in the LEGEND
# instead of showing up as an unexplained grey block. See .fig1a_top_annotation().
.fig1a_na_label <- "not classified"
.fig1a_na_col   <- "#FFFFFF"        # white = absent; never used for a real category

.fig1a_covariate_cols <- function(a) {
  cols <- list()
  # Integrated TCGA class — factor level order is POLE > MMRd > CN-high (serous) > CN-low
  # (endometrioid). Assigned positionally so it tracks the levels the reports build.
  # The `not classified` level (added by .fig1a_top_annotation when values are missing) is
  # held OUT of the positional assignment and always gets white, so a gap can never be
  # mistaken for a class — in TCGA 2013 roughly a THIRD of samples have no subtype, so this
  # is the single largest block of colour in the annotation.
  if ("TCGA_class" %in% names(a) && is.factor(a$TCGA_class)) {
    lv   <- levels(a$TCGA_class)
    real <- setdiff(lv, .fig1a_na_label)
    cc   <- stats::setNames(attend_pal(length(real)), real)     # positional fallback
    # ...then override by NAME wherever the level is one we recognise. The integrated-class
    # bar sits directly above the MMR and Aneuploidy bars, so a reader reads the three as a
    # stack: giving MMRd the MMR bar's orange, CN-high the aneuploidy-high red and CN-low
    # the aneuploidy-low light blue makes the stack agree with itself instead of assigning
    # a fourth unrelated colour to the same concept one row down.
    # Named, not positional, because a positional palette silently reshuffles every colour
    # if attend_tcga_levels ever gains a tier (e.g. POLE, present in TCGA but not ATTEND).
    # POLE takes the highlight-reserved blue (rule [C]) and that is safe here for two
    # independent reasons: ATTEND has no POLE cases, so the level never materialises; and
    # this is a heatmap ANNOTATION BAR, which never carries the per-patient highlight
    # points that rule [C] protects. Kept in the map so adding POLE upstream cannot
    # silently reshuffle the positional fallback.
    known <- c(`POLE`                             = unname(ATTEND_PALETTE["blue"]),
               `MMRd`                             = unname(attend_mmr_cols[["MMRd"]]),
               `Copy-number high (serous-like)`   = attend_aneu_high,
               `Copy-number low (endometrioid)`   = attend_aneu_low)
    hit <- intersect(names(known), real)
    if (length(hit)) cc[hit] <- known[hit]
    if (.fig1a_na_label %in% lv) cc[[.fig1a_na_label]] <- .fig1a_na_col
    cols$TCGA_class <- cc
  }
  # MMR IHC and TP53 — read from the semantic palettes rather than re-spelled here, for the
  # same reason the aneuploidy bars below do: a hex written twice is a hex that drifts. Both
  # take the palette's reference grey (rule [A]), NOT a pale tint — these are
  # real measurements and must not read as blank cells.
  if ("MMR" %in% names(a))  cols$MMR  <- c(MMRp = unname(attend_mmr_cols[["MMRp"]]),
                                           MMRd = unname(attend_mmr_cols[["MMRd"]]))
  if ("TP53" %in% names(a)) cols$TP53 <- c(wt  = unname(attend_tp53_cols[["wt"]]),
                                           mut = unname(attend_tp53_cols[["mut"]]))
  # A PUBLISHED partition drawn beside the recomputed one (report 10: TCGA's own
  # CNA_CLUSTER_K4). It takes the SAME .fig1a_cluster_cols the bottom cluster bar uses, so
  # "recomputed 3" and "published 3" are the same colour and agreement is legible by eye
  # instead of requiring the reader to hold a legend mapping in their head. Left to
  # ComplexHeatmap's defaults this bar would get an unrelated palette and the comparison —
  # the entire point of the panel — would be the hardest thing on the figure to read.
  for (nm in intersect(c("Published_CN_cluster", "CN_cluster_published"), names(a))) {
    lv <- if (is.factor(a[[nm]])) levels(a[[nm]]) else sort(unique(stats::na.omit(as.character(a[[nm]]))))
    if (length(lv)) cols[[nm]] <- .fig1a_cluster_cols[lv]
  }
  # Binary aneuploidy — light blue = low, red = high. Read from attend_aneu_cols rather
  # than re-spelled here, so the heatmap bar and the ggplot boxplots cannot drift apart.
  if ("Aneuploidy_hl" %in% names(a))
    cols$Aneuploidy_hl <- c(`aneu-low`  = unname(attend_aneu_cols[["aneu-low"]]),
                            `aneu-high` = unname(attend_aneu_cols[["aneu-high"]]))
  # Continuous aneuploidy score — DIVERGING light blue -> white -> red, pinned at the SAME
  # high/low split the binary bar uses (.aneu_split_threshold), so white sits exactly on the
  # cutoff and the two aneuploidy bars agree cell by cell. Falls back to a plain two-stop
  # ramp when the threshold is not strictly inside the observed range.
  if ("Aneuploidy" %in% names(a) && is.numeric(a$Aneuploidy)) {
    rng <- range(a$Aneuploidy, na.rm = TRUE)
    lo  <- unname(attend_aneu_cols[["aneu-low"]]); hi <- unname(attend_aneu_cols[["aneu-high"]])
    if (all(is.finite(rng)) && diff(rng) > 0) {
      thr <- .aneu_split_threshold(a$Aneuploidy)
      cols$Aneuploidy <- if (is.finite(thr) && thr > rng[1] && thr < rng[2])
        circlize::colorRamp2(c(rng[1], thr, rng[2]), c(lo, "#F7F7F7", hi))
      else circlize::colorRamp2(rng, c(lo, hi))
    }
  }
  cols
}

# Build the top HeatmapAnnotation for the Fig-1a covariate bars — shared by fig1a_heatmap()
# and fig1a_heatmap_ksweep() so all TCGA reports get the same treatment. Applies three
# cross-report conventions in ONE place:
#   (1) DROP the redundant CN_cluster bar — the bottom Cluster bar already encodes it;
#   (2) DERIVE a binary Aneuploidy_hl (high/low) from the continuous Aneuploidy score. The
#       split uses the pipeline cutoff attend_thresholds$aneuploidy when the score is on the
#       [0,1] master scale, else the cohort median — so it is valid whether the bar came from
#       ASCETS (arm-count scale) or the master aneuploidy score;
#   (3) colour every bar from the readable .fig1a_covariate_cols() palette.
# Returns a column HeatmapAnnotation, or NULL when nothing is left to draw. `ids` = the
# heatmap column order (barcodes). Assumes ComplexHeatmap/circlize (callers guard first).
.fig1a_top_annotation <- function(ann_col, ids) {
  if (is.null(ann_col)) return(NULL)
  a <- ann_col[ids, , drop = FALSE]
  a$CN_cluster <- NULL                                  # (1) redundant with the bottom bar
  # (2) binary aneuploidy high/low from the continuous score (same split as the boxplots)
  if ("Aneuploidy" %in% names(a) && is.numeric(a$Aneuploidy) && any(is.finite(a$Aneuploidy))) {
    thr <- .aneu_split_threshold(a$Aneuploidy)
    a$Aneuploidy_hl <- factor(ifelse(is.na(a$Aneuploidy) | !is.finite(thr), NA_character_,
                                     ifelse(a$Aneuploidy >= thr, "aneu-high", "aneu-low")),
                              levels = c("aneu-low", "aneu-high"))
  }
  # (2b) Make MISSING an explicit, LEGENDED category on TCGA_class rather than an
  # unlabelled grey. Without this, ComplexHeatmap paints NA in its default grey (#BEBEBE),
  # which was almost the same colour as MMRp (#BBBBBB) and aneu-low (#D9D9D9) in the rows
  # directly beneath — so "no subtype", "MMR-proficient" and "aneuploidy-low" all looked
  # alike. In TCGA 2013, 131 of 363 heatmap columns (36%) have no integrated subtype,
  # because it was only assigned to the 232-sample core set, so this is the single biggest
  # block of colour in the annotation and must be self-explanatory.
  if ("TCGA_class" %in% names(a) && any(is.na(a$TCGA_class))) {
    a$TCGA_class <- factor(ifelse(is.na(a$TCGA_class), .fig1a_na_label,
                                  as.character(a$TCGA_class)),
                           levels = c(levels(a$TCGA_class), .fig1a_na_label))
  }
  # bar order: integrated class, MMR, TP53, then binary + continuous aneuploidy last
  ord <- c("TCGA_class", "MMR", "TP53", "Aneuploidy_hl", "Aneuploidy")
  a   <- a[, c(intersect(ord, names(a)), setdiff(names(a), ord)), drop = FALSE]
  a   <- a[, vapply(a, function(x) !all(is.na(x)), logical(1)), drop = FALSE]  # drop all-NA
  if (!ncol(a)) return(NULL)
  # na_col = white so any REMAINING missing value (MMR, aneuploidy) reads as an empty cell
  # rather than as a category; border = TRUE keeps those white cells visible as cells.
  ComplexHeatmap::HeatmapAnnotation(df = a, which = "column",
                                    col = .fig1a_covariate_cols(a),
                                    na_col = .fig1a_na_col, border = TRUE,
                                    annotation_name_side = "left")
}

# ============================================================================
# Cluster labels ordered by SCNA burden — the paper's x-axis order.
# ============================================================================
# WHY THIS EXISTS. In Kandoth et al. Fig. 1a the four cluster blocks run left to right from
# copy-number QUIET to copy-number HIGH: measured on TCGA's own published CNA_CLUSTER_K4 over
# their hg19 segments, mean |log2| per cluster is 0.0043, 0.052, 0.099, 0.216 — monotonic, a
# 50x range, and their "Copy-number high (Serous-like)" subtype is 60/60 inside cluster 4.
# TCGA numbered their clusters by burden, so the figure reads as a gradient.
#
# cutree() does not. It numbers clusters by order of first appearance in the dendrogram, which
# is an artefact of merge order and carries no meaning. Splitting a heatmap on those labels
# puts the blocks in an arbitrary sequence, which is why a faithful-looking reproduction still
# does not look like the paper: same data, same clusters, wrong order.
#
# Relabelling by ascending mean |value| fixes the axis AND buys two things: cluster k is the
# copy-number-high group by construction, and cluster numbers become comparable ACROSS cohorts
# (ATTEND cluster 4 and TCGA cluster 4 both mean "most altered"), which a cutree label never is.
#
# Burden is computed on the matrix passed in. Callers pass the DISPLAYED (genome-wide,
# continuous) matrix, because that is the quantity the paper's own ordering tracks and it is
# the matrix both the single-k panel and the k-sweep already hold.
#
# Returns the relabelled vector AND the old->new map, because a caller that has already named
# a cluster by some other rule — cnv_high_cluster_burden() on the peak matrix, say — must
# translate that name rather than assume it is now k. Base R, so it is unit-testable without
# the tidyverse stack.
relabel_clusters_by_burden <- function(clusters, mat) {
  old <- as.character(clusters)
  nm  <- names(clusters)
  lv  <- unique(stats::na.omit(old))
  if (!length(lv) || is.null(mat) || is.null(rownames(mat)) || is.null(nm))
    return(list(cluster = clusters, map = stats::setNames(lv, lv), burden = NULL))

  # mean |value| per cluster, over the samples the matrix actually carries. A cluster with no
  # row in `mat` gets NA and is ordered last, so it can never silently take slot 1.
  b <- vapply(lv, function(l) {
    ids <- nm[!is.na(old) & old == l]
    ids <- intersect(ids, rownames(mat))
    if (!length(ids)) return(NA_real_)
    mean(abs(as.matrix(mat[ids, , drop = FALSE])), na.rm = TRUE)
  }, numeric(1))

  ord <- order(b, na.last = TRUE)
  map <- stats::setNames(as.character(seq_along(lv)), lv[ord])
  new <- factor(unname(map[old]), levels = as.character(seq_along(lv)))
  names(new) <- nm
  list(cluster = new, map = map, burden = stats::setNames(unname(b), lv))
}

# Render one tumour x chromosomal-location SCNA heatmap in the TCGA Fig. 1a
# layout with ComplexHeatmap: rows = features split by chromosome (chr1..X down
# the left edge), columns = tumours split by copy-number cluster with the cluster
# bar beneath, red=amp/blue=del. `ann_col` (optional) is a data.frame keyed by
# the SAME rownames as `mat` (barcodes) carrying the existing 09/09b annotation
# bars (TCGA_class, MMR, TP53, Aneuploidy); drawn as a top annotation. Knit-safe:
# with ComplexHeatmap/circlize absent it draws the cluster dendrogram instead.
fig1a_heatmap <- function(mat, feature_pos, clusters, ann_col = NULL,
                          value_type = c("thresholded", "continuous"), main = NULL) {
  value_type <- match.arg(value_type)
  if (is.null(mat) || !nrow(mat) || !ncol(mat)) { message("fig1a_heatmap(): empty matrix"); return(invisible(NULL)) }
  # Drop samples with no cluster assignment so column_split never receives NA
  keep <- rownames(mat) %in% names(clusters)
  if (!all(keep)) {
    message("fig1a_heatmap(): dropping ", sum(!keep), " sample(s) with no cluster assignment")
    mat <- mat[keep, , drop = FALSE]
  }
  if (!nrow(mat)) { message("fig1a_heatmap(): no samples with a cluster"); return(invisible(NULL)) }
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    message("fig1a_heatmap(): ComplexHeatmap/circlize not installed — install for the faithful Fig-1a; ",
            "drawing a plain dendrogram fallback.")
    hc <- tryCatch(stats::hclust(stats::dist(mat)), error = function(e) NULL)
    if (!is.null(hc)) plot(hc, labels = FALSE, main = main %||% "CNV clustering (fallback)")
    return(invisible(NULL))
  }

  # rows = features present in mat, ordered by chromosomal location
  fp  <- feature_pos[feature_pos$feature %in% colnames(mat), , drop = FALSE]
  if (!nrow(fp)) { message("fig1a_heatmap(): no feature positions resolved"); return(invisible(NULL)) }
  M   <- t(mat[, fp$feature, drop = FALSE])            # features (rows) x samples (cols)
  row_chr <- droplevels(fp$chrom)

  # columns split by cluster (paper groups tumours into the k clusters)
  ids     <- colnames(M)
  cl_vec  <- factor(unname(clusters[ids]))
  cols    <- .fig1a_cluster_cols[levels(cl_vec)]; names(cols) <- levels(cl_vec)

  # color scale: continuous red-white-blue, or discrete GISTIC {-2..+2}
  if (value_type == "continuous") {
    lim    <- stats::quantile(abs(M), 0.98, na.rm = TRUE); lim <- if (is.finite(lim) && lim > 0) lim else 1
    col_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("blue", "white", "red"))
    leg_title <- "SCNA (log2)"
  } else {
    col_fun <- circlize::colorRamp2(c(-2, -1, 0, 1, 2),
                                    c("navy", "royalblue", "white", "firebrick2", "darkred"))
    leg_title <- "GISTIC call"
  }

  # bottom = cluster bar (TCGA colors); top = the ATTEND covariate bars, built centrally by
  # .fig1a_top_annotation() (drops the redundant CN_cluster, adds a binary aneuploidy bar,
  # and applies the readable covariate palette — same treatment for every TCGA report).
  bottom <- ComplexHeatmap::HeatmapAnnotation(
    Cluster = cl_vec, col = list(Cluster = cols),
    annotation_name_side = "left", which = "column")
  top <- .fig1a_top_annotation(ann_col, ids)

  ht <- ComplexHeatmap::Heatmap(
    M, name = leg_title, col = col_fun,
    row_split = row_chr, cluster_row_slices = FALSE, cluster_rows = FALSE,
    row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 7),
    column_split = cl_vec, cluster_column_slices = FALSE, cluster_columns = TRUE,
    show_row_names = FALSE, show_column_names = FALSE,
    top_annotation = top, bottom_annotation = bottom,
    column_title = main %||% "SCNAs by tumour (columns) x chromosomal location (rows) — TCGA Fig. 1a",
    row_title_side = "left", border = TRUE,
    heatmap_legend_param = list(title = leg_title))
  ComplexHeatmap::draw(ht, merge_legend = TRUE)
  invisible(ht)
}

# Fig-1a chromosome heatmap for a WHOLE k-sweep in ONE figure. The SCNA matrix is drawn
# once (features x tumours, rows split by chromosome, columns ordered by the supplied
# dendrogram `hc` — here the Euclidean/Ward one from cluster_arm_matrix), and beneath it a
# STACK of cluster bars, one per k in `k_range`: bar "k=NN" colours each tumour by
# cutree(hc, NN). Reading down the bars shows how the same tumours partition as k grows
# (2 -> 8), all against the same chromosomal SCNA background. `ann_col` optional top
# annotation (TCGA_class/MMR/TP53/Aneuploidy), keyed by the rownames of `mat` as in
# fig1a_heatmap(). Because the columns follow the ACTUAL dendrogram (not
# ComplexHeatmap's default column clustering), the bars are contiguous blocks. Knit-safe:
# ComplexHeatmap/circlize absent -> the dendrogram is drawn with the k-cut rectangles.
fig1a_heatmap_ksweep <- function(mat, feature_pos, hc, k_range = 2:8, ann_col = NULL,
                                 value_type = c("thresholded", "continuous"), main = NULL) {
  value_type <- match.arg(value_type)
  if (is.null(mat) || !nrow(mat) || !ncol(mat)) { message("fig1a_heatmap_ksweep(): empty matrix"); return(invisible(NULL)) }
  if (is.null(hc)) { message("fig1a_heatmap_ksweep(): no dendrogram"); return(invisible(NULL)) }

  # keep only samples shared by the matrix and the dendrogram; clamp k to [2, n-1]
  ids <- intersect(rownames(mat), hc$labels)
  if (length(ids) < 3) { message("fig1a_heatmap_ksweep(): < 3 shared samples"); return(invisible(NULL)) }
  ks  <- k_range[k_range >= 2 & k_range <= (length(ids) - 1)]
  if (!length(ks)) { message("fig1a_heatmap_ksweep(): no valid k in range"); return(invisible(NULL)) }

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    message("fig1a_heatmap_ksweep(): ComplexHeatmap/circlize not installed — drawing the ",
            "clustering dendrogram with k-cut rectangles as a fallback.")
    plot(hc, labels = FALSE, main = main %||% "CNV clustering (Ward.D2) — k-sweep",
         xlab = "", sub = "")
    for (k in ks) stats::rect.hclust(hc, k = k, border = .fig1a_cluster_cols[as.character(k)])
    return(invisible(NULL))
  }

  # rows = features present in mat, ordered by chromosomal location
  fp <- feature_pos[feature_pos$feature %in% colnames(mat), , drop = FALSE]
  if (!nrow(fp)) { message("fig1a_heatmap_ksweep(): no feature positions resolved"); return(invisible(NULL)) }
  M       <- t(mat[ids, fp$feature, drop = FALSE])   # features (rows) x samples (cols)
  row_chr <- droplevels(fp$chrom)

  # color scale: continuous red-white-blue, or discrete GISTIC {-2..+2}
  if (value_type == "continuous") {
    lim     <- stats::quantile(abs(M), 0.98, na.rm = TRUE); lim <- if (is.finite(lim) && lim > 0) lim else 1
    col_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("blue", "white", "red"))
    leg_title <- "SCNA (log2)"
  } else {
    col_fun <- circlize::colorRamp2(c(-2, -1, 0, 1, 2),
                                    c("navy", "royalblue", "white", "firebrick2", "darkred"))
    leg_title <- "GISTIC call"
  }

  # bottom = one cluster bar per k (cutree(hc, k)); colours from the shared TCGA palette.
  # Every k is relabelled by ascending SCNA burden, exactly as the single-k panel is, so the
  # bars read as a gradient and the stack is internally consistent: raw cutree labels would
  # make "cluster 1" mean a different thing on every row of the sweep.
  bars <- as.data.frame(lapply(ks, function(k) {
    cut_k <- stats::cutree(hc, k = k)[ids]
    relabel_clusters_by_burden(cut_k, mat[ids, , drop = FALSE])$cluster
  }))
  names(bars) <- paste0("k=", ks); rownames(bars) <- ids
  col_list <- stats::setNames(lapply(names(bars), function(nm) {
    lv <- levels(bars[[nm]]); stats::setNames(.fig1a_cluster_cols[lv], lv)
  }), names(bars))
  bottom <- ComplexHeatmap::HeatmapAnnotation(df = bars, col = col_list, which = "column",
                                              annotation_name_side = "left",
                                              gap = grid::unit(0.4, "mm"),
                                              simple_anno_size = grid::unit(3, "mm"))

  # top = ATTEND covariate bars, built centrally by .fig1a_top_annotation() (drops the
  # redundant CN_cluster, adds the binary aneuploidy bar, readable palette) — same as fig1a_heatmap()
  top <- .fig1a_top_annotation(ann_col, ids)

  ht <- ComplexHeatmap::Heatmap(
    M, name = leg_title, col = col_fun,
    row_split = row_chr, cluster_row_slices = FALSE, cluster_rows = FALSE,
    row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 7),
    cluster_columns = stats::as.dendrogram(hc), column_dend_reorder = FALSE,
    show_row_names = FALSE, show_column_names = FALSE,
    top_annotation = top, bottom_annotation = bottom,
    column_title = main %||% paste0("SCNAs by tumour (columns, dendrogram order) x chromosomal ",
                                    "location (rows); cluster bars k = ", min(ks), "..", max(ks)),
    row_title_side = "left", border = TRUE,
    heatmap_legend_param = list(title = leg_title))
  ComplexHeatmap::draw(ht, merge_legend = TRUE)
  invisible(ht)
}

# ============================================================================
# Manual patient highlight — shared by aneu_point_layers() (report 08) AND every
# other per-patient point/box plot. Centralising the id
# matching here means the reports cannot drift on which id spaces are tried or on
# the highlight colours. Base-R sourceable; ggplot only touched inside highlight_points().
# ============================================================================

# Normalise an id vector to a canonical string form so a patient matches regardless of how
# a loader typed its id column. Excel/readxl reads a numeric subject id as a DOUBLE, and the
# round-trip can surface as "785", "785.0", or " 785 " depending on the path — none of which
# string-equal each other. This trims whitespace and strips a trailing ".0" (a double rendered
# as text), leaving alphanumeric ids ("21S188") untouched. Applied to BOTH sides of the match.
.norm_id <- function(x) sub("\\.0+$", "", trimws(as.character(x)))

# --- Resolving a configured highlight id into every id space -----------------------
# THE BUG THIS EXISTS FOR. The highlight ids are configured in whatever space the analyst
# knows a sample by: `21S188` is a sequencing BARCODE. Most figures, though, are built from
# the master, which is keyed on `pid` and carries no barcode column (collapse_pid() prefixes
# every other dataset's columns, so even a surviving id becomes `aneu__ID`, which is not one
# of the names highlight_group_of() looks for). A literal-string match therefore found
# nothing and the overlay silently did nothing on the majority of plots — the worst failure
# mode an overlay has, because the figure still renders, just without the mark.
#
# So the ids are translated through the SAME crosswalks the join uses. attend_plots.R cannot
# build them (they need tidyverse, and this file is base-R sourceable so the unit tests run
# without the stack), so the crosswalks are pushed IN: build_crosswalks() registers
# barcode->pid and image->pid on the fresh path, and get_master() re-registers them from the
# written intermediates on the cached path. Every report picks this up with no call-site
# change, which is the same reason get_master() exists at all.
#
# Kept in an environment rather than a bare global so register_*() cannot be defeated by
# `<<-` resolving into the wrong frame when a report sources this file from inside a function.
.attend_hl_env <- new.env(parent = emptyenv())
.attend_hl_env$xwalk <- data.frame(id = character(0), pid = character(0),
                                   stringsAsFactors = FALSE)

#' Register id -> pid maps for highlight resolution.
#'
#' Accepts any number of maps, each either a two-column data.frame (the column literally
#' named `pid` is the pid side; the other is the foreign id — this is the shape
#' build_barcode_pid()/build_image_pid() return) or a named character vector (names = foreign
#' id, values = pid). Replaces the stored table by default so a rebuild cannot accumulate
#' stale ids; `append = TRUE` to add a map from somewhere else.
#' @return invisibly, the number of id -> pid pairs now registered.
register_highlight_xwalk <- function(..., append = FALSE) {
  rows <- list()
  for (m in list(...)) {
    if (is.null(m) || !length(m)) next
    if (is.data.frame(m)) {
      if (ncol(m) < 2L) next
      pid_col <- if ("pid" %in% names(m)) "pid" else names(m)[[2L]]
      id_col  <- setdiff(names(m)[seq_len(2L)], pid_col)
      id_col  <- if (length(id_col)) id_col[[1L]] else names(m)[[1L]]
      rows[[length(rows) + 1L]] <- data.frame(id  = .norm_id(m[[id_col]]),
                                              pid = .norm_id(m[[pid_col]]),
                                              stringsAsFactors = FALSE)
    } else if (!is.null(names(m))) {
      rows[[length(rows) + 1L]] <- data.frame(id  = .norm_id(names(m)),
                                              pid = .norm_id(unname(m)),
                                              stringsAsFactors = FALSE)
    }
  }
  new <- if (length(rows)) do.call(rbind, rows) else
    data.frame(id = character(0), pid = character(0), stringsAsFactors = FALSE)
  keep <- !is.na(new$id) & !is.na(new$pid) & nzchar(new$id) & nzchar(new$pid)
  new  <- new[keep, , drop = FALSE]
  if (append) new <- rbind(.attend_hl_env$xwalk, new)
  .attend_hl_env$xwalk <- unique(new)
  invisible(nrow(.attend_hl_env$xwalk))
}

#' Expand configured ids to every id that names the SAME patient.
#'
#' One hop each way is enough: both crosswalks point AT `pid`, so
#' barcode -> pid -> {every barcode, every image id} reaches the whole equivalence class.
#' Returns the configured ids unchanged when nothing is registered, so the pre-registration
#' behaviour (literal match only) is preserved and this is never a new failure mode.
highlight_expand_ids <- function(ids, xwalk = .attend_hl_env$xwalk) {
  ids <- .norm_id(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids) || is.null(xwalk) || !nrow(xwalk)) return(unique(ids))
  pids <- unique(c(ids[ids %in% xwalk$pid], xwalk$pid[xwalk$id %in% ids]))
  unique(c(ids, pids, xwalk$id[xwalk$pid %in% pids]))
}

# Return the highlight GROUP NAME for each row of `df` (NA where none). Matches against
# whichever id column df carries — pid / ID / barcode / Tumor_Sample_Barcode / image_id /
# patient_id — after expanding the configured ids through the registered crosswalks
# (highlight_expand_ids), so a barcode-configured sample is still found on a pid-keyed
# master and vice-versa. The FIRST group in highlight$groups wins, so a patient listed in
# two groups keeps the earlier one; with a single group configured that is a no-op, but it
# is the invariant a second group would depend on, so it stays.
# Matching is TYPE-AGNOSTIC: normalised-string equality OR numeric equality (so a numeric id
# column matches string ids and vice-versa). Knit-safe: no id column or no config -> all-NA.
highlight_group_of <- function(df, id_cols = c("pid", "ID", "barcode", "Tumor_Sample_Barcode",
                                               "image_id", "patient_id"),
                               highlight = if (exists("attend_highlight")) attend_highlight else NULL) {
  df   <- as.data.frame(df)
  out  <- rep(NA_character_, nrow(df))
  cand <- intersect(id_cols, names(df))
  if (is.null(highlight) || !length(highlight$groups) || !length(cand) || !nrow(df)) return(out)
  for (nm in names(highlight$groups)) {
    ids <- highlight$groups[[nm]]$ids
    if (!length(ids)) next
    # Crosswalk expansion, not a bare .norm_id(): see highlight_expand_ids() above for why
    # a literal match silently missed 21S188 on every master-derived figure.
    ids_norm <- highlight_expand_ids(ids)
    ids_num  <- suppressWarnings(as.numeric(ids_norm))
    hit <- Reduce(`|`, lapply(cand, function(cc) {
      col_norm <- .norm_id(df[[cc]])
      m <- col_norm %in% ids_norm                                  # normalised-string match
      col_num <- suppressWarnings(as.numeric(col_norm))            # numeric fallback (both sides parse)
      m | (!is.na(col_num) & col_num %in% ids_num[!is.na(ids_num)])
    }))
    out[hit & is.na(out)] <- nm
  }
  out
}

# One-line coverage diagnostic — how many of each highlight group's ids resolved in `df`, and
# which did NOT. Print it beside a plot (`cat(highlight_coverage(df))`) to see, on real data,
# whether the ids match the id column at all. Pure text; safe to leave in a report.
highlight_coverage <- function(df, highlight = if (exists("attend_highlight")) attend_highlight else NULL,
                               id_cols = c("pid", "ID", "barcode", "Tumor_Sample_Barcode",
                                           "image_id", "patient_id")) {
  df   <- as.data.frame(df)
  cand <- intersect(id_cols, names(df))
  if (is.null(highlight) || !length(highlight$groups))
    return("highlight: no groups configured\n")
  if (!length(cand)) return(paste0("highlight: no id column found (looked for ",
                                    paste(id_cols, collapse = "/"), ")\n"))
  present <- unique(unlist(lapply(cand, function(cc) .norm_id(df[[cc]]))))
  present_num <- suppressWarnings(as.numeric(present))
  # Scored per CONFIGURED id but resolved through the crosswalks, so the line reads
  # "polipo: 1/1 matched" on a pid-keyed master even though 21S188 is a barcode. Reporting
  # the raw configured id as missing while the expansion matched would be a false alarm.
  lines <- vapply(names(highlight$groups), function(nm) {
    cfg  <- highlight$groups[[nm]]$ids
    hit  <- vapply(cfg, function(one) {
      ids <- highlight_expand_ids(one)
      any(ids %in% present |
          (!is.na(suppressWarnings(as.numeric(ids))) &
           suppressWarnings(as.numeric(ids)) %in% present_num[!is.na(present_num)]))
    }, logical(1))
    miss <- cfg[!hit]
    sprintf("  %s: %d/%d matched%s", nm, sum(hit), length(cfg),
            if (length(miss)) paste0(" | missing: ", paste(miss, collapse = ", ")) else "")
  }, character(1))
  paste0("highlight coverage (id cols: ", paste(cand, collapse = "/"), ")\n",
         paste(lines, collapse = "\n"), "\n")
}

# Named colour vector c(<group> = "#RRGGBB", …) from the highlight config, in the
# declared group order. Empty when nothing is configured.
highlight_palette <- function(highlight = if (exists("attend_highlight")) attend_highlight else NULL) {
  if (is.null(highlight) || !length(highlight$groups)) return(character(0))
  vapply(highlight$groups, function(g) g$color %||% "#111111", character(1))
}

# ggplot layers that draw a per-patient plot's points with the manual highlight applied:
# highlighted patients are drawn LARGER, ON TOP, in their group colour with a labelled
# "highlight" legend; everyone else in the base layer. Points are PARTITIONED (a highlighted
# patient appears ONLY in the overlay), so no point is ever drawn twice — the invariant
# aneu_point_layers() relies on, which matters once points are jittered (two independent jitters
# would double-place a patient).
#
# The highlight rides its OWN colour scale so it gets a legend. When the base already colours by a
# categorical (`base_colour`), the two colour scales are separated with ggnewscale::new_scale_colour()
# — so highlight_points must OWN the base colour scale (pass `base_palette` for a manual one, else a
# default hue scale is used) to control the ggnewscale ordering; the caller must NOT add its own base
# colour scale. When the base has no colour (`base_colour = NULL`), the highlight simply claims the
# single free colour scale — no ggnewscale needed. If ggnewscale is unavailable but the base owns
# colour, it degrades to a constant-colour overlay (visible, no legend) rather than breaking.
#
# `x`/`y` are the plot's axis COLUMN NAMES. `base_palette` = named colours for `base_colour`'s levels
# (NULL -> default hue); `base_name`/`base_na` title & NA colour of the base legend. `legend_title`
# names the highlight legend. `jitter_width = 0` -> geom_point (scatters); > 0 -> geom_jitter (fixed
# seed). Knit-safe: ggplot2 absent, or no rows, -> empty list.
# Usage (free colour): `p + highlight_points(df, "group", "aneuploidy_value")`
# Usage (coloured base): `p + highlight_points(df, "x", "y", base_colour = "cls", base_palette = pal)`
highlight_points <- function(df, x, y, base_colour = NULL, base_palette = NULL,
                             base_name = NULL, base_na = "grey70",
                             jitter_width = 0.15, seed = 1,
                             base_size = 0.7, base_alpha = 0.6, base_const = "grey25",
                             hl_size = 3, hl_stroke = 1, hl_alpha = 1, legend_title = "highlight",
                             highlight = if (exists("attend_highlight")) attend_highlight else NULL,
                             id_cols = c("pid", "ID", "barcode", "Tumor_Sample_Barcode",
                                         "image_id", "patient_id")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(list())
  df <- as.data.frame(df)
  if (!nrow(df)) return(list())
  df$.hl <- highlight_group_of(df, id_cols, highlight)
  # attend_highlight_show = FALSE (attend_classes.R) suppresses the MARK, not the patient.
  # This has to clear .hl BEFORE the split below: `base` is the is.na(.hl) rows, so gating
  # any later would drop the highlighted patient out of the figure altogether — and at every
  # call site that passes points = FALSE to attend_box(), highlight_points() IS the points
  # layer, so that silently deletes a patient from the plotted cohort.
  # exists() because attend_plots.R is sourced in a bootstrap env with no attend_classes.R.
  if (exists("attend_highlight_show") && !isTRUE(attend_highlight_show))
    df$.hl <- NA_character_
  base   <- df[is.na(df$.hl), , drop = FALSE]
  hl     <- df[!is.na(df$.hl), , drop = FALSE]
  pal    <- highlight_palette(highlight)
  bname  <- if (is.null(base_name)) ggplot2::waiver() else base_name

  # one geom constructor for every layer: jittered (fixed seed) or plain points
  mk <- if (jitter_width > 0)
    function(...) ggplot2::geom_jitter(
      position = ggplot2::position_jitter(width = jitter_width, height = 0, seed = seed), ...)
  else function(...) ggplot2::geom_point(...)

  # --- base layer (+ its own colour scale when base_colour is set, so ggnewscale ordering is ours)
  if (!is.null(base_colour)) {
    layers <- list(mk(data = base,
                      mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = .data[[base_colour]]),
                      size = base_size, alpha = base_alpha))
    layers <- c(layers, list(
      if (!is.null(base_palette)) ggplot2::scale_colour_manual(values = base_palette, na.value = base_na, name = bname)
      else                        ggplot2::scale_colour_hue(na.value = base_na, name = bname)))
  } else {
    layers <- list(mk(data = base, mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]]),
                      size = base_size, alpha = base_alpha, colour = base_const))
  }
  if (!nrow(hl)) return(layers)
  hl$.hl <- factor(hl$.hl, levels = names(pal))

  has_ngs <- requireNamespace("ggnewscale", quietly = TRUE)
  if (!is.null(base_colour) && !has_ngs) {
    # base owns the colour scale and we cannot open a second one: constant-colour overlay (no legend)
    for (nm in names(pal)) {
      sub <- hl[hl$.hl == nm, , drop = FALSE]; if (!nrow(sub)) next
      layers <- c(layers, list(mk(data = sub, mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]]),
                                  colour = pal[[nm]], size = hl_size, stroke = hl_stroke, alpha = hl_alpha)))
    }
    return(layers)
  }
  if (!is.null(base_colour)) layers <- c(layers, list(ggnewscale::new_scale_colour()))

  # highlight overlay on its own colour scale -> labelled legend
  c(layers, list(
    mk(data = hl, mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = .data[[".hl"]]),
       size = hl_size, stroke = hl_stroke, alpha = hl_alpha),
    ggplot2::scale_colour_manual(values = pal, name = legend_title, drop = FALSE, na.translate = FALSE)))
}

# Point layers for the aneuploidy-by-cluster / -by-class boxplots (report 08). Encodes two
# things on the jittered points WITHOUT touching the boxplot's cluster/class fill:
#   (1) SHAPE marks aneuploidy-HIGH tumours (triangle) vs low (circle), split by
#       .aneu_split_threshold() so it agrees with the heatmap's binary aneuploidy bar;
#   (2) a MANUAL highlight subset (attend_highlight$groups — patient `pid`s or barcodes, matched
#       against whichever id column the plot carries) is redrawn ON TOP in the group's colour,
#       larger, so a hand-picked cohort stands out. Colour is reserved for the highlight so it
#       never competes with the aneu-high shape.
# Returns a list of ggplot2 layers/scales to ADD to an existing ggplot (`p + aneu_point_layers(df)`).
# Empty highlight -> just the high/low shape split. `df` must carry `value_col` and, for the
# highlight, an id column (pid/ID/barcode). Base-R friendly; ggplot2 attached by the reports.
aneu_point_layers <- function(df, value_col = "aneuploidy_value",
                              highlight = if (exists("attend_highlight")) attend_highlight else NULL,
                              seed = 1) {
  df <- as.data.frame(df)
  v   <- suppressWarnings(as.numeric(df[[value_col]]))
  thr <- .aneu_split_threshold(v)
  df$.aneu_high <- factor(ifelse(!is.na(v) & is.finite(thr) & v >= thr, "aneu-high", "aneu-low"),
                          levels = c("aneu-low", "aneu-high"))
  # manual highlight: match ids against any plausible id column present in df
  # (shared with the other reports via highlight_group_of()/highlight_palette()).
  df$.hl  <- highlight_group_of(df, highlight = highlight)
  col_map <- highlight_palette(highlight)
  pj    <- ggplot2::position_jitter(width = 0.15, height = 0, seed = seed)
  base  <- df[is.na(df$.hl), , drop = FALSE]
  hlpts <- df[!is.na(df$.hl), , drop = FALSE]
  layers <- list(
    ggplot2::geom_point(data = base, ggplot2::aes(shape = .data[[".aneu_high"]]),
                        position = pj, size = 1.7, alpha = 0.65, colour = "grey20"),
    ggplot2::scale_shape_manual(values = c(`aneu-low` = 16, `aneu-high` = 17),
                                name = "aneuploidy", drop = FALSE)
  )
  if (nrow(hlpts)) {
    hlpts$.hl <- factor(hlpts$.hl, levels = names(col_map))
    layers <- c(layers, list(
      ggplot2::geom_point(data = hlpts,
                          ggplot2::aes(shape = .data[[".aneu_high"]], colour = .data[[".hl"]]),
                          position = pj, size = 3, stroke = 1),
      ggplot2::scale_colour_manual(values = col_map, name = "highlight",
                                   drop = FALSE, na.translate = FALSE)
    ))
  }
  layers
}

# TCGA Fig. 1b — Kaplan-Meier of progression-free survival for each copy-number
# cluster. Maps clustered barcodes -> pid (via the report-01 crosswalk) -> the
# clinical PFS columns, recodes the free-text event with recode_event(), fits
# survfit(Surv ~ cnv_cluster), draws the curves (ggsurvfit if available, else
# base), and returns the log-rank p (survdiff). Knit-safe: NULL if survival is
# absent, the PFS columns are missing, or < 2 clusters carry an event.
km_by_cluster <- function(assignment, cw, master,
                          time_col  = attend_cols$surv_time,
                          event_col = attend_cols$surv_event) {
  if (is.null(assignment) || is.null(cw) || is.null(master)) return(invisible(NULL))
  if (!requireNamespace("survival", quietly = TRUE)) { message("km_by_cluster(): survival not installed"); return(invisible(NULL)) }
  if (!all(c(time_col, event_col) %in% names(master))) { message("km_by_cluster(): PFS columns absent"); return(invisible(NULL)) }

  d <- merge(assignment, cw, by.x = "ID", by.y = "barcode")
  d <- merge(d, master[, c("pid", time_col, event_col)], by = "pid")
  d$time  <- suppressWarnings(as.numeric(d[[time_col]]))
  d$event <- recode_event(d[[event_col]])                 # free text -> 0/1
  d$cnv_cluster <- factor(d$cnv_cluster)
  d <- d[!is.na(d$time) & !is.na(d$event) & !is.na(d$cnv_cluster), , drop = FALSE]
  if (nrow(d) < 2 || nlevels(droplevels(d$cnv_cluster)) < 2) { message("km_by_cluster(): too few groups"); return(invisible(NULL)) }
  d$cnv_cluster <- droplevels(d$cnv_cluster)

  fit  <- survival::survfit(survival::Surv(time, event) ~ cnv_cluster, data = d)
  sd   <- survival::survdiff(survival::Surv(time, event) ~ cnv_cluster, data = d)
  pval <- stats::pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)

  if (requireNamespace("ggsurvfit", quietly = TRUE)) {
    print(ggsurvfit::ggsurvfit(fit) +
            ggplot2::labs(y = "Progression-free survival", x = "Survival (months)",
                          title = sprintf("PFS by copy-number cluster (log-rank P = %.4g)", pval)))
  } else {
    plot(fit, col = .fig1a_cluster_cols[levels(d$cnv_cluster)], lwd = 2,
         xlab = "Survival (months)", ylab = "Progression-free survival",
         main = sprintf("PFS by copy-number cluster (log-rank P = %.4g)", pval))
    legend("bottomleft", legend = paste("Cluster", levels(d$cnv_cluster)),
           col = .fig1a_cluster_cols[levels(d$cnv_cluster)], lwd = 2, bty = "n")
  }
  invisible(list(p_logrank = pval, n = nrow(d), fit = fit))
}

# --- report 06: recurrent SCNA frequency plots --------------------------------

#' Genomic midpoint of each peak, for the x axis.
.peak_pos <- function(peaks) {
  data.frame(peak_id = peaks$peak_id,
             chrom   = factor(peaks$chrom,
                              levels = c(as.character(1:22), "X", "Y")),
             pos     = (peaks$wide_start + peaks$wide_end) / 2,
             direction = peaks$direction, stringsAsFactors = FALSE)
}

#' GISTIC-style mirrored frequency plot: amplifications up, deletions down.
scna_mirror_plot <- function(freq_tbl, peaks, title = NULL) {
  if (is.null(freq_tbl) || !nrow(freq_tbl)) {
    message("scna_mirror_plot(): no data — skipping."); return(NULL)
  }
  d <- merge(freq_tbl, .peak_pos(peaks), by = "peak_id")
  if (!nrow(d)) { message("scna_mirror_plot(): no peaks matched."); return(NULL) }
  d$signed <- ifelse(d$direction == "del", -d$freq, d$freq)

  ggplot2::ggplot(d, ggplot2::aes(x = pos, y = signed, fill = direction)) +
    ggplot2::geom_col(width = 2e6) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::facet_grid(scna_group ~ chrom, scales = "free_x", space = "free_x") +
    # amp/del is a DIVERGING direction encoding, not one of the semantic categorical
    # palettes, so it may use the highlight-reserved blue: a frequency panel has no
    # per-patient points for a highlight to disappear into, and a diverging pair needs
    # opposing hues. red/blue about zero, matching the SCNA heatmap body.
    ggplot2::scale_fill_manual(values = c(amp = unname(ATTEND_PALETTE["red"]),
                                          del = unname(ATTEND_PALETTE["blue"]))) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(abs(x) * 100, "%")) +
    ggplot2::labs(x = "Genomic position", y = "Altered (del below / amp above)",
                  title = title) +
    attend_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.spacing.x = ggplot2::unit(0.05, "lines"))
}

#' Frequency difference (aneuploidy-high minus -low) along the genome, per MMR stratum.
scna_delta_plot <- function(freq_tbl, peaks, title = NULL) {
  if (is.null(freq_tbl) || !nrow(freq_tbl)) {
    message("scna_delta_plot(): no data — skipping."); return(NULL)
  }
  d <- freq_tbl
  d$mmr <- sub("-.*$", "", as.character(d$scna_group))
  d$hl  <- sub("^.*-", "", as.character(d$scna_group))
  w <- stats::reshape(d[, c("peak_id", "mmr", "hl", "freq")],
                      idvar = c("peak_id", "mmr"), timevar = "hl",
                      direction = "wide")
  if (!all(c("freq.high", "freq.low") %in% names(w))) {
    message("scna_delta_plot(): need both high and low groups — skipping."); return(NULL)
  }
  w$delta <- w$freq.high - w$freq.low
  d2 <- merge(w, .peak_pos(peaks), by = "peak_id")
  if (!nrow(d2)) { message("scna_delta_plot(): no peaks matched."); return(NULL) }

  ggplot2::ggplot(d2, ggplot2::aes(x = pos, y = delta, colour = mmr)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
    ggplot2::scale_colour_manual(values = attend_mmr_cols) +
    ggplot2::labs(x = "Genomic position",
                  y = "Frequency difference (high - low)",
                  colour = "MMR", title = title) +
    attend_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.spacing.x = ggplot2::unit(0.05, "lines"))
}

#' Composite serous-like panel score by group — report 06's primary figure.
# `ylab` is an argument because there are now TWO panel scores on this plot type and
# they measure different things: the pre-specified TCGA panel (Family A) and the
# data-driven one (Family B). A hardcoded "TCGA UCEC panel" label on the Family B
# figure would assert exactly the provenance that panel does not have.
panel_score_plot <- function(scores, group, title = NULL,
                             ylab = "Fraction of TCGA UCEC panel altered (0-1)") {
  if (is.null(scores) || !length(scores)) {
    message("panel_score_plot(): no scores — skipping."); return(NULL)
  }
  d <- data.frame(score = as.numeric(scores), grp = group, stringsAsFactors = FALSE)
  d <- d[!is.na(d$grp) & !is.na(d$score), , drop = FALSE]
  if (!nrow(d)) { message("panel_score_plot(): all NA — skipping."); return(NULL) }

  # No `fill = grp`: the group is already the x position, so colouring by it spends the
  # palette on information the axis carries — which is why the old version had to switch
  # its own legend off. A single neutral fill, and attend_box() picks box vs median-crossbar
  # from each group's n (these are SCNA strata; MMRd-high can be a handful of patients).
  ggplot2::ggplot(d, ggplot2::aes(x = grp, y = score)) +
    attend_box(d, x = "grp", y = "score", point_size = 1.2, point_alpha = 0.75) +
    ggplot2::labs(x = NULL, y = ylab, title = title) +
    attend_theme() +
    ggplot2::theme(legend.position = "none")
}
