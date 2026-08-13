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

# TCGA Fig. 1a color for the copy-number clusters (green/purple/blue/red …),
# matching the paper's Cluster bar. Covers up to 8 clusters so the k = 2..8 sweep
# (fig1a_heatmap_ksweep()) always has a distinct colour per cluster; the standard
# pipeline fixes k at 4 (Fig. 2b). k > 8 would need more colors.
.fig1a_cluster_cols <- c("1" = "#2CA02C", "2" = "#7B4FA3", "3" = "#3B6FB6", "4" = "#D62728",
                         "5" = "#FF7F0E", "6" = "#8C564B", "7" = "#17BECF", "8" = "#BCBD22")

# Project-wide qualitative palette — the brewer-free replacement for
# scale_*_brewer(palette = "Set1") / RColorBrewer::brewer.pal(). Reports use
#   scale_colour_manual(values = attend_pal(n), aesthetics = c("colour", "fill"))
# instead, so no plot depends on RColorBrewer being loadable.
# Colours are the Okabe-Ito colour-blind-safe set, the same one .fig1a_covariate_cols()
# already uses for the heatmap annotation bars, so categorical colours agree across reports.
# `n` > 8 recycles with a warning — a categorical plot needing more than 8 colours should
# use a different encoding anyway. Pass `names` to get a named vector for a manual scale.
attend_pal <- function(n = 8, names = NULL) {
  okabe_ito <- c("#0072B2",  # blue
                 "#D55E00",  # vermillion
                 "#009E73",  # bluish green
                 "#CC79A7",  # reddish purple
                 "#E69F00",  # orange
                 "#56B4E9",  # sky blue
                 "#F0E442",  # yellow
                 "#999999")  # grey
  if (!is.null(names)) n <- length(names)
  if (n > length(okabe_ito))
    warning("attend_pal(): ", n, " colours requested but only ", length(okabe_ito),
            " distinct ones exist — recycling.", call. = FALSE)
  out <- rep(okabe_ito, length.out = max(1, n))
  if (!is.null(names)) stats::setNames(out, names) else out
}

# --- Semantic colour palettes (single source of truth) --------------------
# Promoted from .fig1a_covariate_cols() so the ggplot reports and the Fig-1a
# heatmap use the SAME hex per concept. Keyed by every label spelling in use
# across reports so a single manual scale covers heatmap bars and ggplots alike.
attend_mmr_cols <- c(MMRp = "#999999", MMRd = "#E69F00",
                     `MMR proficient` = "#999999", `MMR deficient` = "#E69F00")
# Aneuploidy — LOW light blue, HIGH red. Chosen to read the same way as the SCNA
# heatmap body itself (blue = quiet/loss end, red = altered/gain end) so the annotation
# bar and the matrix beneath it do not tell opposite colour stories. Deliberately NOT
# #B2182B (that is TP53-abnormal, which sits in the adjacent Fig-1a bar) and not the
# Okabe-Ito orange (MMRd).
attend_aneu_low  <- "#A6CEE3"   # light blue
attend_aneu_high <- "#E31A1C"   # red
attend_aneu_cols <- c(`aneuploidy-low` = attend_aneu_low, `aneuploidy-high` = attend_aneu_high,
                      `aneu-low` = attend_aneu_low, `aneu-high` = attend_aneu_high,
                      `MMRd aneuploidy-low` = attend_aneu_low,
                      `MMRd aneuploidy-high` = attend_aneu_high)
attend_tp53_cols <- c(`TP53-normal` = "#AAAAAA", `TP53-abnormal` = "#B2182B",
                      wt = "#AAAAAA", mut = "#B2182B")

# Treatment response (add_response_class): responder bluish-green, non-responder reddish
# purple. Both are Okabe-Ito, and both are deliberately outside every palette above --
# response is plotted on the same cell-type panels as aneuploidy, so reusing the aneuploidy
# blue/red would make two different groupings look like the same one at a glance.
attend_resp_cols <- c(`responder` = "#009E73", `non-responder` = "#CC79A7",
                      `MMRd responder` = "#009E73", `MMRd non-responder` = "#CC79A7")

# Project-wide ggplot theme. One base_size (11) and legend/caption conventions so
# fonts and spacing stop drifting between reports (was theme_minimal(11|12) and
# theme_bw()). Call sites still add theme(legend.position = ...) on top as needed.
attend_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(size = base_size - 1),
      strip.text    = ggplot2::element_text(face = "bold"),
      legend.title  = ggplot2::element_text(size = base_size - 1),
      panel.grid.minor = ggplot2::element_blank())
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
    pal  <- c("#CC79A7", "#009E73", "#D55E00", "#0072B2")  # purple, green, vermillion, blue
    cc   <- stats::setNames(pal[seq_along(real)], real)
    if (.fig1a_na_label %in% lv) cc[[.fig1a_na_label]] <- .fig1a_na_col
    cols$TCGA_class <- cc
  }
  # MMR IHC — proficient MID grey, deficient orange. Deliberately #999999 (the Okabe-Ito
  # grey), NOT a pale tint: it is a real measurement and must not read as a blank cell.
  if ("MMR" %in% names(a))  cols$MMR  <- c(MMRp = "#999999", MMRd = "#E69F00")
  # TP53 — wild-type mid grey (same reasoning as MMRp), mutant deep red
  if ("TP53" %in% names(a)) cols$TP53 <- c(wt = "#AAAAAA", mut = "#B2182B")
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
# dendrogram `hc` — here the Pearson/Ward one from cluster_arm_matrix), and beneath it a
# STACK of cluster bars, one per k in `k_range`: bar "k=NN" colours each tumour by
# cutree(hc, NN). Reading down the bars shows how the same tumours partition as k grows
# (2 -> 8), all against the same chromosomal SCNA background. `ann_col` optional top
# annotation (TCGA_class/MMR/TP53/Aneuploidy), keyed by the rownames of `mat` as in
# fig1a_heatmap(). Because the columns follow the ACTUAL Pearson/Ward dendrogram (not
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
            "Pearson/Ward dendrogram with k-cut rectangles as a fallback.")
    plot(hc, labels = FALSE, main = main %||% "CNV clustering (Pearson/Ward) — k-sweep",
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

  # bottom = one cluster bar per k (cutree(hc, k)); colours from the shared TCGA palette
  bars <- as.data.frame(lapply(ks, function(k) factor(unname(stats::cutree(hc, k = k)[ids]))))
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
    column_title = main %||% paste0("SCNAs by tumour (columns, Pearson/Ward order) x chromosomal ",
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

# Return the highlight GROUP NAME for each row of `df` (NA where none). Matches the configured
# ids against whichever id column df carries — pid / ID / barcode / image_id / patient_id (the
# imaging space is included because the highlight ids may be imaging-side). The FIRST group in
# highlight$groups wins, so a patient listed in two groups keeps the earlier one — this is why
# attend_highlight lists `polipo` before `cohort`, so 21S188 can never be pulled into the cohort.
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
    ids_norm <- .norm_id(ids)
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
  lines <- vapply(names(highlight$groups), function(nm) {
    ids  <- .norm_id(highlight$groups[[nm]]$ids)
    hit  <- ids %in% present |
            (!is.na(suppressWarnings(as.numeric(ids))) &
             suppressWarnings(as.numeric(ids)) %in% present_num[!is.na(present_num)])
    miss <- highlight$groups[[nm]]$ids[!hit]
    sprintf("  %s: %d/%d matched%s", nm, sum(hit), length(ids),
            if (length(miss)) paste0(" | missing: ", paste(miss, collapse = ", ")) else "")
  }, character(1))
  paste0("highlight coverage (id cols: ", paste(cand, collapse = "/"), ")\n",
         paste(lines, collapse = "\n"), "\n")
}

# Named colour vector c(<group> = "#RRGGBB", …) from the highlight config, in the
# declared group order (polipo, cohort). Empty when nothing is configured.
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
    ggplot2::scale_fill_manual(values = c(amp = "#D55E00", del = "#0072B2")) +
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
panel_score_plot <- function(scores, group, title = NULL) {
  if (is.null(scores) || !length(scores)) {
    message("panel_score_plot(): no scores — skipping."); return(NULL)
  }
  d <- data.frame(score = as.numeric(scores), grp = group, stringsAsFactors = FALSE)
  d <- d[!is.na(d$grp) & !is.na(d$score), , drop = FALSE]
  if (!nrow(d)) { message("panel_score_plot(): all NA — skipping."); return(NULL) }

  ggplot2::ggplot(d, ggplot2::aes(x = grp, y = score, fill = grp)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 1.8, alpha = 0.85) +
    ggplot2::labs(x = NULL, y = "Fraction of TCGA UCEC panel altered",
                  title = title) +
    attend_theme() +
    ggplot2::theme(legend.position = "none")
}
