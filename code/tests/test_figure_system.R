# Base-R unit tests (no tidyverse) — runnable in the bootstrap env.
source(file.path("code", "attend_plots.R"))  # base-R-sourceable: no top-level library()

## --- semantic palettes exist and carry the promoted Fig-1a hexes ---
stopifnot(is.character(attend_mmr_cols),
          attend_mmr_cols[["MMRp"]] == "#999999",
          attend_mmr_cols[["MMRd"]] == "#E69F00",
          attend_mmr_cols[["MMR proficient"]] == "#999999",
          attend_mmr_cols[["MMR deficient"]] == "#E69F00")

stopifnot(identical(attend_aneu_low, "#A6CEE3"),   # light blue = aneuploidy-low
          identical(attend_aneu_high, "#E31A1C"),  # red        = aneuploidy-high
          attend_aneu_high != attend_tp53_cols[["mut"]],  # must not clash with the adjacent
          attend_aneu_high != attend_mmr_cols[["MMRd"]],  # Fig-1a annotation bars
          is.character(attend_aneu_cols),
          attend_aneu_cols[["aneuploidy-low"]]  == attend_aneu_low,
          attend_aneu_cols[["aneuploidy-high"]] == attend_aneu_high,
          attend_aneu_cols[["aneu-low"]]  == attend_aneu_low,
          attend_aneu_cols[["aneu-high"]] == attend_aneu_high)

stopifnot(is.character(attend_tp53_cols),
          attend_tp53_cols[["TP53-normal"]]   == "#AAAAAA",
          attend_tp53_cols[["TP53-abnormal"]] == "#B2182B")

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

## --- attend_aneu_cols carries the MMRd-subclass spellings (Phase 1.5) ---
stopifnot(attend_aneu_cols[["MMRd aneuploidy-low"]]  == attend_aneu_low,
          attend_aneu_cols[["MMRd aneuploidy-high"]] == attend_aneu_high)

cat("test_figure_system: ALL PASS\n")
