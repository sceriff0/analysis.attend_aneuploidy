#!/usr/bin/env Rscript
# Test the GISTIC -maxseg hypothesis: how many samples in the combined .seg
# exceed the per-sample segment cap and are therefore dropped by GISTIC.
#   Rscript code/diagnose_maxseg.R
suppressPackageStartupMessages({ library(tidyverse); library(here); library(data.table) })
source(here("code", "attend_classes.R"))

out <- here(attend_cnv$seg$gistic_seg_out)
stopifnot("combined .seg missing — run write_gistic_seg() first" = file.exists(out))
comb <- fread(out, header = TRUE) |> as_tibble()

per <- comb |> count(Sample, name = "n_seg") |> arrange(desc(n_seg))
cat(sprintf("samples total                : %d\n", nrow(per)))
for (cap in c(2000, 2500, 5000, 10000, 25000, 50000)) {
  cat(sprintf("samples with <= %6d segments : %d   (GISTIC keeps these at -maxseg %d)\n",
              cap, sum(per$n_seg <= cap), cap))
}
cat(sprintf("\nsegment-count distribution: min %d / median %.0f / mean %.0f / max %d\n",
            min(per$n_seg), median(per$n_seg), mean(per$n_seg), max(per$n_seg)))
cat("\nTop 5 most-segmented samples:\n"); print(utils::head(per, 5))
cat("\nIf '<= 2000' is ~25, that confirms -maxseg 2000 is dropping the rest.\n")
