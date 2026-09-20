#!/usr/bin/env Rscript
# =============================================================================
# preflight_gistic.R — why the GISTIC scripts fail while every report knits fine
#
# WHY THIS EXISTS. The reports and the GISTIC scripts read DISJOINT inputs. Reports consume
# GISTIC's OUTPUT (data/gistic/, data/ascets/, the written master); the prep scripts consume
# its INPUT (data/seg/, output/gistic_input/). No report reads data/seg/ at all — verified by
# grep, not assumed — so the entire GISTIC chain can be broken while all twelve reports knit
# to completion. "The .Rmd work, the .R don't" is therefore the EXPECTED symptom of a missing
# data/seg/, not a puzzle about Rscript.
#
# Each prep step also stops at its FIRST failure, so a run tells you one missing thing per
# attempt. This checks everything and reports all of it at once, with the fix beside each.
#
# RUN (repo root, on the cluster):   Rscript code/preflight_gistic.R
# Exits 0 when the pipeline can run, 1 when a REQUIRED item is missing.
# =============================================================================
suppressPackageStartupMessages({
  ok_pkg <- vapply(c("tidyverse", "here", "data.table", "fs"),
                   function(p) requireNamespace(p, quietly = TRUE), logical(1))
})
if (!all(ok_pkg)) {
  cat("FAIL  packages: missing ", paste(names(ok_pkg)[!ok_pkg], collapse = ", "), "\n",
      "      fix: renv::restore(), or install them into the R the cluster job uses.\n", sep = "")
  quit(status = 1)
}
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(data.table); library(fs)
})

fails <- 0L
say <- function(state, label, detail = "", fix = "") {
  if (identical(state, "FAIL")) fails <<- fails + 1L
  cat(sprintf("%-5s %-26s %s\n", state, label, detail))
  if (nzchar(fix)) cat("      fix: ", fix, "\n", sep = "")
}

cat("GISTIC preflight — repo root resolved by here(): ", here(), "\n",
    "R ", R.version$major, ".", R.version$minor,
    " | wd ", getwd(), "\n\n", sep = "")

# ⚠️ IS here() EVEN THE REPO? here() locates the project by walking UP for a .Rproj / .here /
# DESCRIPTION marker. Started from somewhere outside the project — a home directory, or an
# interactive R session that never cd'd — it finds no marker and, rather than failing, FALLS
# BACK TO THE WORKING DIRECTORY and returns a confidently wrong root. Every source() below
# then resolves cleanly to a path that does not exist, and the run emits seven identical
# "cannot open the connection" FAILs that read like a code problem and are not one. Checking
# for one file the repo must contain turns that into the one sentence it actually is.
if (!file.exists(file.path(here(), "code", "attend_classes.R"))) {
  cat("FAIL  not the repo root        here() = ", here(), "\n",
      "      code/attend_classes.R is not there, so this is not the ATTEND checkout.\n",
      "      here() falls back to the working directory when it finds no project marker,\n",
      "      so every path below would resolve to a file that does not exist.\n",
      "      fix: cd into the repo and run with Rscript, not by pasting into R:\n",
      "             cd ~/workflowR/attend_aneuploidy\n",
      "             Rscript code/preflight_gistic.R\n", sep = "")
  quit(status = 1)
}

# Sourcing is itself a check: an error here is a code problem, not a data one, and it is the
# one failure mode that WOULD also break the reports.
srcs <- c("attend_classes.R", "attend_scna.R", "attend_harmonise.R", "attend_io.R",
          "attend_plots.R", "load_wes_results.R", "load_clinical.R")
for (f in srcs) {
  e <- tryCatch({ source(here("code", f)); NULL }, error = function(e) conditionMessage(e))
  if (is.null(e)) say("ok", paste0("source ", f)) else
    say("FAIL", paste0("source ", f), e, "this breaks the reports too — fix the code first.")
}
if (fails) { cat("\nStopping: the helper layer does not load.\n"); quit(status = 1) }

cat("\n-- INPUTS THE SCRIPTS NEED AND THE REPORTS DO NOT --------------------------\n")

# 1. data/seg/ — the raw per-sample segments. THE usual answer.
seg_dir <- here("data", attend_cnv$seg$dir)
segs    <- if (dir.exists(seg_dir))
  list.files(seg_dir, pattern = "\\.seg$", recursive = TRUE, full.names = TRUE) else character(0)
# Braces required. At top level an `if (x) a` with `else` on the NEXT line is an R SYNTAX
# error — R closes the if-statement at the newline and then meets a stray `else`. It parses
# fine inside {} or a function body, which is why it survives a quick visual check and dies
# only when the file is actually parsed. Same trap report 11 documents for its dip label.
if (!dir.exists(seg_dir)) {
  say("FAIL", "data/seg/", paste0("no folder at ", seg_dir),
      "sync the DRAGEN per-sample .seg exports here. NO REPORT READS THIS, which is why the site still knits without it.")
} else if (!length(segs)) {
  say("FAIL", "data/seg/", paste0("folder exists, 0 *.seg inside (", seg_dir, ")"),
      "the exports did not land — check the transfer, and that the files end in .seg")
} else {
  say("ok", "data/seg/", paste0(length(segs), " .seg file(s)"))
}

# 2. the master — needed for the GROUP and LOO segs, not for the pooled one.
master <- tryCatch(read_intermediate("attend_master_joined"), error = function(e) NULL)
if (is.null(master) || !nrow(master)) {
  say("FAIL", "master table", "absent or empty",
      "knit report 01, or Rscript code/build_master.R")
} else {
  say("ok", "master table", paste0(nrow(master), " patients x ", ncol(master), " cols"))
}

# 3. the barcode -> pid crosswalk.
cw <- tryCatch(build_barcode_pid(load_gianlu_clinical_data(), make_id_cfg()),
               error = function(e) NULL)
if (is.null(cw) || !nrow(cw)) {
  say("FAIL", "barcode->pid crosswalk", "empty",
      "check load_gianlu_clinical_data() and make_id_cfg()")
} else {
  say("ok", "barcode->pid crosswalk", paste0(nrow(cw), " pairs, ",
                                             dplyr::n_distinct(cw$pid), " patients"))
}

# 4. GROUP SIZES among the barcodes the .seg files actually carry. This is what decides
#    whether each per-group GISTIC run is even possible, and it is far cheaper to read here
#    than after a multi-hour run that dies.
if (length(segs) && !is.null(master) && !is.null(cw)) {
  g <- tryCatch({
    m   <- add_scna_group(add_molecular_classes(master))
    bc  <- seg_barcodes()
    tab <- table(factor(as.character(barcode_scna_group(bc, m, cw)),
                        levels = attend_scna$group_levels), useNA = "ifany")
    list(n_bc = length(bc), tab = tab)
  }, error = function(e) NULL)
  if (is.null(g)) {
    say("FAIL", "scna_group over segs", "could not be computed",
        "usually a column still on its attend_cols placeholder")
  } else {
    say("ok", "scna_group over segs", paste0(g$n_bc, " seg barcodes classified"))
    cat("\n"); print(g$tab); cat("\n")
    thin <- names(g$tab)[g$tab > 0 & g$tab < 3 & !is.na(names(g$tab))]
    if (length(thin))
      say("warn", "thin group(s)", paste(thin, collapse = ", "),
          "GISTIC needs a handful of samples; these runs may start and die, leaving a folder holding only gistic_inputs.mat")
    zero <- names(g$tab)[g$tab == 0 & !is.na(names(g$tab))]
    if (length(zero))
      say("FAIL", "empty group(s)", paste(zero, collapse = ", "),
          "no .seg carries these samples — their attend_<group>.seg cannot be written, so run_one() will skip them and NO FOLDER will appear")
  }
}

cat("\n-- THE SEGS THE PREP STEP WRITES ------------------------------------------\n")
gi <- here("output", "gistic_input")
if (!dir.exists(gi)) {
  say("warn", "output/gistic_input/", "absent — prep will create it",
      "MODE=groups Rscript code/prep_gistic_group_segs.R")
} else {
  want <- c("attend_all_segments.seg", paste0("attend_", c("mmrp_high", "mmrp_low",
                                                           "mmrd_low", "mmrd_high"), ".seg"))
  got  <- file.exists(file.path(gi, want))
  for (i in seq_along(want))
    say(if (got[i]) "ok" else "warn", want[i],
        if (got[i]) "present" else "missing",
        if (got[i]) "" else "run_one() SKIPS a group whose .seg is missing and creates no folder — this is why data/gistic/<group>/ can be absent entirely")
  n_loo <- length(list.files(gi, pattern = "^attend_mmrd_high_drop_.*\\.seg$"))
  say(if (n_loo > 0) "ok" else "warn", "LOO segs", paste0(n_loo, " found"),
      if (n_loo > 0) "" else "MODE=loo Rscript code/prep_gistic_group_segs.R")
}

cat("\n-- GISTIC RUNS ALREADY ON DISK --------------------------------------------\n")
gx <- here("data", attend_cnv$gistic$dir)
rd <- setNames(file.path(gx, c("all", "mmrp_high", "mmrp_low", "mmrd_low", "mmrd_high")),
               c("all", "mmrp_high", "mmrp_low", "mmrd_low", "mmrd_high"))
print(gistic_run_inventory(rd), row.names = FALSE)

# The peak-gene lists are what the report-07 clustering restricts on. A run can be "complete"
# by file count and still yield zero parsed gene symbols, which silently un-restricts the
# clustering — so parse them here rather than discovering it in a knit.
apk <- tryCatch(.gistic_peak_genes(find_gistic_files()$amp_genes), error = function(e) character(0))
dpk <- tryCatch(.gistic_peak_genes(find_gistic_files()$del_genes), error = function(e) character(0))
cat("\n")
if (length(union(apk, dpk)) == 0) {
  say("FAIL", "peak-gene parse", "0 symbols from amp_genes + del_genes",
      "report 07 then clusters ALL genes instead of the significant peaks. Send the first 6 lines of the amp_genes file to check the parser.")
} else {
  say("ok", "peak-gene parse", paste0(length(apk), " amp + ", length(dpk),
                                      " del = ", length(union(apk, dpk)), " unique peak genes"))
}

cat("\n-- GISTIC ITSELF (what run_gistic.sh needs) --------------------------------\n")
sif <- Sys.getenv("GISTIC_SIF"); bin <- Sys.getenv("GISTIC_BIN", "gistic2")
mod <- Sys.getenv("MODULE");     ref <- Sys.getenv("REFGENE",
  "refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat")
if (nzchar(sif)) {
  say(if (file.exists(sif)) "ok" else "FAIL", "GISTIC_SIF", sif,
      if (file.exists(sif)) "" else "path does not exist")
} else if (nzchar(mod)) {
  say("ok", "MODULE", mod, "")
} else {
  found <- nzchar(Sys.which(bin))
  say(if (found) "ok" else "FAIL", "GISTIC_BIN", paste0(bin, if (found) paste0(" -> ", Sys.which(bin)) else " not on PATH"),
      if (found) "" else "set GISTIC_SIF=<image>.sif, or MODULE=gistic/2.0.23, or put gistic2 on PATH")
}
say(if (file.exists(ref)) "ok" else "warn", "REFGENE", ref,
    if (file.exists(ref)) "" else "relative to the repo root; override with REFGENE=<path>")

cat("\n", strrep("-", 76), "\n", sep = "")
if (fails) {
  cat(fails, " REQUIRED item(s) missing — fix those, then:\n",
      "  MODE=groups Rscript code/prep_gistic_group_segs.R\n",
      "  MODE=groups bash code/run_gistic.sh\n",
      "  MODE=loo    Rscript code/prep_gistic_group_segs.R\n",
      "  MODE=loo    bash code/run_gistic.sh\n", sep = "")
  quit(status = 1)
}
cat("All required inputs present. Run:\n",
    "  MODE=groups Rscript code/prep_gistic_group_segs.R && MODE=groups bash code/run_gistic.sh\n",
    "  MODE=loo    Rscript code/prep_gistic_group_segs.R && MODE=loo    bash code/run_gistic.sh\n", sep = "")
