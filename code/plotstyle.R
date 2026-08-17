# Canonical publication figure style (R / ggplot2).
#
# Source of truth: ~/.claude/plotstyle/plotstyle.R
# Copied into repos as plotstyle.R; do not edit the copy, edit the source.
# Deliberately depends on ggplot2 + base grid only, so it runs anywhere.
#
# Usage
#   source("plotstyle.R")
#   p <- ggplot(df, aes(x, y, colour = group)) + geom_point() +
#        scale_colour_pub(df$group, figdir = "figures") + theme_pub()
#   pub_save(p, "figures/fig2_expression", width = "single", data = df)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# ---------------------------------------------------------------------------
# Palette -- Okabe & Ito (2008), colourblind-safe, greyscale-legible.
# Order is locked and must match plotstyle.py exactly.
# ---------------------------------------------------------------------------
OKABE_ITO <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
               "#56B4E9", "#D55E00", "#F0E442", "#000000")
GREY <- "#4D4D4D"; LIGHT_GREY <- "#BFBFBF"
SEQUENTIAL <- "viridis"; DIVERGING <- "RdBu"

# Journal column widths, millimetres. Must match plotstyle.py.
WIDTHS_MM <- c(single = 85, onehalf = 114, double = 180)
MAX_HEIGHT_MM <- 230
BASE_PT <- 7; TICK_PT <- 6; PANEL_PT <- 8

# ggplot linewidths are in mm; convert from points.
lw <- function(pt) pt / .pt

theme_pub <- function(base_size = BASE_PT, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      line = element_line(linewidth = lw(0.5), colour = "black"),
      rect = element_rect(fill = "white", colour = NA),
      text = element_text(colour = "black", size = base_size),
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      axis.line        = element_line(linewidth = lw(0.5), colour = "black"),
      axis.ticks       = element_line(linewidth = lw(0.5), colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      axis.text        = element_text(size = TICK_PT, colour = "black"),
      axis.title       = element_text(size = base_size),
      plot.title       = element_text(size = base_size, hjust = 0,
                                      margin = margin(b = 3)),
      plot.subtitle    = element_text(size = TICK_PT, hjust = 0, colour = GREY),
      plot.tag         = element_text(size = PANEL_PT, face = "bold"),
      plot.tag.position = c(0, 1),
      legend.key       = element_blank(),
      legend.background = element_blank(),
      legend.title     = element_text(size = base_size),
      legend.text      = element_text(size = TICK_PT),
      legend.key.size  = unit(3, "mm"),
      legend.margin    = margin(0, 0, 0, 0),
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size, hjust = 0,
                                      margin = margin(b = 2)),
      plot.margin      = margin(2, 2, 2, 2)
    )
}

# ---------------------------------------------------------------------------
# Locked category -> colour mapping, shared with the Python module through the
# same .palette.json file. This is what stops a group changing colour between
# figure 1 and figure 3.
# ---------------------------------------------------------------------------
lock_palette <- function(categories, figdir = "figures", name = "default") {
  categories <- as.character(unique(categories))
  path <- file.path(figdir, ".palette.json")
  store <- list()
  if (file.exists(path)) {
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      store <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    } else {
      warning("jsonlite not installed: palette cannot be shared with Python plots")
    }
  }
  mapping <- if (!is.null(store[[name]])) unlist(store[[name]]) else character(0)
  for (cat in categories) {
    if (!cat %in% names(mapping)) {
      if (length(mapping) >= length(OKABE_ITO))
        stop(sprintf("%d categories exceeds the 8-colour palette; facet, group, or encode with shape",
                     length(mapping) + 1L))
      mapping[cat] <- OKABE_ITO[length(mapping) + 1L]
    }
  }
  store[[name]] <- as.list(mapping)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
    writeLines(jsonlite::toJSON(store, auto_unbox = TRUE, pretty = TRUE), path)
  }
  mapping[categories]
}

scale_colour_pub <- function(categories, figdir = "figures", name = "default", ...)
  scale_colour_manual(values = lock_palette(categories, figdir, name), ...)
scale_color_pub <- scale_colour_pub
scale_fill_pub <- function(categories, figdir = "figures", name = "default", ...)
  scale_fill_manual(values = lock_palette(categories, figdir, name), ...)

# ---------------------------------------------------------------------------
# Composition. Uses patchwork when installed, otherwise a base-grid layout so
# panels stay the same size and tags stay put.
# ---------------------------------------------------------------------------
pub_compose <- function(plots, ncol = NULL, nrow = NULL, tags = TRUE) {
  n <- length(plots)
  if (is.null(ncol) && is.null(nrow)) ncol <- ceiling(sqrt(n))
  if (is.null(ncol)) ncol <- ceiling(n / nrow)
  if (is.null(nrow)) nrow <- ceiling(n / ncol)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    pw <- patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow)
    if (tags && n > 1) pw <- pw + patchwork::plot_annotation(tag_levels = "a")
    return(pw)
  }
  if (tags && n > 1)
    plots <- Map(function(p, t) p + labs(tag = t), plots, letters[seq_len(n)])
  # Convert to gtables and force every panel to the same geometry. Without
  # this, a panel carrying a legend or a longer axis label gets a smaller
  # plotting region than its neighbours -- the ragged multi-panel look.
  grobs <- lapply(plots, ggplot2::ggplotGrob)
  nw <- unique(vapply(grobs, function(g) length(g$widths), integer(1)))
  nh <- unique(vapply(grobs, function(g) length(g$heights), integer(1)))
  if (length(nw) == 1L && length(nh) == 1L) {
    maxw <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
    maxh <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$heights))
    grobs <- lapply(grobs, function(g) { g$widths <- maxw; g$heights <- maxh; g })
  } else {
    warning("panels have incompatible gtable layouts; sizes may differ")
  }
  function() {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow, ncol)))
    for (i in seq_len(n)) {
      pushViewport(viewport(layout.pos.row = ceiling(i / ncol),
                            layout.pos.col = (i - 1L) %% ncol + 1L))
      grid.draw(grobs[[i]])
      popViewport()
    }
    popViewport()
  }
}

# ---------------------------------------------------------------------------
# Output: PDF + PNG + SVG + source-data CSV, from one call.
# ---------------------------------------------------------------------------
pub_save <- function(plot, path, width = "single", height_mm = NULL,
                     aspect = 0.75, nrow = 1, ncol = 1, data = NULL,
                     formats = c("pdf", "png", "svg"), dpi = 600) {
  w_mm <- if (is.character(width)) unname(WIDTHS_MM[width]) else as.numeric(width)
  if (is.na(w_mm)) stop("width must be single/onehalf/double or a number in mm")
  if (is.null(height_mm)) height_mm <- (w_mm / ncol) * aspect * nrow
  if (height_mm > MAX_HEIGHT_MM)
    stop(sprintf("height %.0fmm exceeds one page (%.0fmm); split into a supplementary figure",
                 height_mm, MAX_HEIGHT_MM))
  stem <- sub("\\.(pdf|png|svg)$", "", path)
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  win <- w_mm / 25.4; hin <- height_mm / 25.4
  draw <- function() if (is.function(plot)) plot() else print(plot)
  written <- character(0)
  for (fmt in formats) {
    out <- paste0(stem, ".", fmt)
    switch(fmt,
      pdf = if (capabilities("cairo")) grDevices::cairo_pdf(out, width = win, height = hin)
            else grDevices::pdf(out, width = win, height = hin, useDingbats = FALSE),
      png = grDevices::png(out, width = win, height = hin, units = "in",
                           res = dpi, type = if (capabilities("cairo")) "cairo" else "quartz"),
      svg = grDevices::svg(out, width = win, height = hin),
      stop("unsupported format: ", fmt))
    draw(); grDevices::dev.off()
    written <- c(written, out)
  }
  if (!is.null(data)) {
    csv <- paste0(stem, "_source_data.csv")
    utils::write.csv(data, csv, row.names = FALSE)
    written <- c(written, csv)
  }
  invisible(written)
}
