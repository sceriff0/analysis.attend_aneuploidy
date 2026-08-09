#!/usr/bin/env Rscript
# Localize the 200-.seg -> 25-sample loss in the GISTIC-input pipeline.
# Run from the ATTEND repo root, where the .seg live:
#   Rscript code/diagnose_seg_pipeline.R
# Counts distinct samples at each hop: data/seg/*.seg  ->  write_gistic_seg's
# combined .seg  ->  GISTIC's all_thresholded.by_genes.txt, and explains every
# .seg that write_gistic_seg would silently skip (unresolved columns / empty).

suppressPackageStartupMessages({ library(tidyverse); library(here); library(fs); library(data.table) })
source(here("code", "attend_classes.R"))   # attend_cnv config

sc <- attend_cnv$seg
folder <- here("data", sc$dir)

cat("\n================ SEG -> GISTIC PIPELINE ================\n")
cat("seg folder (attend_cnv$seg$dir): ", folder, "\n", sep = "")
cat("exists: ", dir_exists(folder), "\n", sep = "")

# ---- hop 1: raw .seg files on disk ----------------------------------------
files <- if (dir_exists(folder)) dir_ls(folder, recurse = TRUE, type = "file", glob = "*.seg") else character(0)
cat(sprintf("\n[hop 1] .seg files under %s (recursive): %d\n", folder, length(files)))
if (length(files) == 0) {
  cat("  -> No .seg here. Your 200 files must be elsewhere; point attend_cnv$seg$dir\n",
      "     at them, or move/symlink them into this folder.\n", sep = "")
  quit(save = "no")
}

# ---- replicate write_gistic_seg's per-file resolution to see what it skips --
pick <- function(nms, cands) { h <- intersect(cands, nms); if (length(h)) h[[1]] else NA_character_ }
status <- map_dfr(files, function(p) {
  dt <- tryCatch(fread(p, header = TRUE, na.strings = c(".", "", "NA"), fill = TRUE, nrows = 5) |> as_tibble(),
                 error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0)
    return(tibble(file = path_file(p), outcome = "EMPTY/UNREADABLE", missing = NA_character_))
  need <- c(chrom = pick(names(dt), sc$chrom_col), start = pick(names(dt), sc$start_col),
            end   = pick(names(dt), sc$end_col),   value = pick(names(dt), sc$value_col))
  miss <- names(need)[is.na(need)]
  tibble(file = path_file(p),
         outcome = if (length(miss)) "SKIPPED (unresolved cols)" else "kept",
         missing = if (length(miss)) paste(miss, collapse = ",") else NA_character_)
})

cat("\n[hop 2] write_gistic_seg() per-file outcome:\n")
print(count(status, outcome, name = "n_files"))
skipped <- filter(status, outcome != "kept")
if (nrow(skipped) > 0) {
  cat("\n  Example skipped files and which required column could not be matched:\n")
  print(utils::head(skipped, 12))
  cat("\n  -> If these are real .seg, their headers differ from the candidate lists in\n",
      "     attend_cnv$seg (chrom/start/end/value_col). Show me one header and I'll\n",
      "     add the variant. Peek a header with:\n",
      "       readLines('", as.character(files[[which(status$outcome!='kept')[1]]]), "', n = 1)\n", sep = "")
}

# ---- hop 3: the combined GISTIC input, if already written ------------------
out <- here(sc$gistic_seg_out)
cat(sprintf("\n[hop 3] combined GISTIC input %s: exists=%s\n", out, file_exists(out)))
if (file_exists(out)) {
  comb <- fread(out, header = TRUE) |> as_tibble()
  cat(sprintf("  distinct samples in combined .seg: %d\n", dplyr::n_distinct(comb$Sample)))
}

# ---- hop 4: what GISTIC actually emitted -----------------------------------
gfolder <- here("data", attend_cnv$gistic$dir)
gt <- if (dir_exists(gfolder))
  dir_ls(gfolder, recurse = TRUE, type = "file", glob = attend_cnv$gistic$all_thresholded_glob) else character(0)
cat(sprintf("\n[hop 4] all_thresholded.by_genes.txt: %s\n",
            if (length(gt)) as.character(gt[[1]]) else "NOT FOUND"))
if (length(gt)) {
  hdr <- fread(gt[[1]], header = TRUE, nrows = 0)
  cat(sprintf("  sample columns (cols minus 3 meta): %d\n", max(0, ncol(hdr) - 3)))
}
cat("\n=======================================================\n")
cat("Read the hops top-down: the first place the count drops below ~200 is the culprit.\n")
