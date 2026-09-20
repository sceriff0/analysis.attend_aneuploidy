#!/usr/bin/env Rscript
# =============================================================================
# diagnose_gistic_runs.R — why SOME GISTIC group runs produced nothing
#
# WHY THIS EXISTS, AND WHY diagnose_maxseg.R IS NOT ENOUGH. diagnose_maxseg.R answers the
# -maxseg question for the POOLED seg only. But GISTIC is run once per scna_group, and the
# groups are small and unequal — MMRd-high is 9 patients. A cap that leaves the pooled run
# comfortable can leave a 9-patient group with too few samples to analyse, and GISTIC then
# dies with a MATLAB error that never names the cause. This checks every run separately,
# on BOTH sides: the input it was given, and the output it produced.
#
# HOW AN EMPTY OUTPUT FOLDER IS READ. run_gistic.sh's run_one() guards the missing seg
# BEFORE mkdir -p:
#     [ -f "$seg_file" ] || { echo "skip: ... missing"; return 0; }
#     mkdir -p "$out_dir"
# so the two states are distinguishable on disk and mean opposite things:
#   NO DIRECTORY          -> the .seg was missing; the fix is upstream, in
#                            prep_gistic_group_segs.R (MODE=groups).
#   DIRECTORY, NO FILES   -> the .seg was found and GISTIC ITSELF failed. Re-running the
#                            same command changes nothing. This is the maxseg / too-few-
#                            samples path.
#   DIRECTORY, SOME FILES -> GISTIC ingested the seg (sample_seg_counts.txt and
#                            gistic_inputs.mat are written early) and then crashed deeper in.
#   all_lesions.conf_*    -> the run is usable; every consumer keys on this file.
#
# RUN (repo root, on the cluster):  Rscript code/diagnose_gistic_runs.R
# Honours MAXSEG from the environment so it tests the cap run_gistic.sh will actually use.
# Exits 0 when every expected run is complete, 1 otherwise.
# =============================================================================
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tibble); library(readr)
})

# ⚠️ Run from the repo root, with Rscript. here() falls back to the working directory when it
# finds no project marker, so from a home directory it returns a confidently wrong root and
# the source() below dies with "cannot open the connection" — a path problem wearing a code
# problem's error message. Same guard, and the same reasoning, as preflight_gistic.R.
if (!file.exists(file.path(here(), "code", "attend_classes.R")))
  stop("not the ATTEND repo root: here() = ", here(),
       "\n  fix: cd ~/workflowR/attend_aneuploidy && Rscript code/diagnose_gistic_runs.R",
       call. = FALSE)

source(here("code", "attend_classes.R"))

# The cap GISTIC will apply. run_gistic.sh defaults to 46000; taking it from the same env
# var means this script tests the cap you are about to run with, not a hardcoded guess.
MAXSEG <- suppressWarnings(as.numeric(Sys.getenv("MAXSEG", "46000")))

segdir <- dirname(here(attend_cnv$seg$gistic_seg_out))
outdir <- here("data", attend_cnv$gistic$dir)

# The pooled run is named in config (`pooled_dir`); the four strata are the scna_group
# levels put through the same tokeniser run_gistic.sh and prep_gistic_group_segs.R use, so
# this list cannot drift from the folders they actually write.
runs <- tibble(
  run = c(attend_cnv$gistic$pooled_dir, "mmrp_high", "mmrp_low", "mmrd_low", "mmrd_high"),
  seg = c(here(attend_cnv$seg$gistic_seg_out),
          file.path(segdir, paste0("attend_", c("mmrp_high", "mmrp_low",
                                                "mmrd_low", "mmrd_high"), ".seg")))
)

# Per-sample segment counts for one .seg. Only the Sample column is read: these files carry
# one row per segment and can run to millions of rows, and nothing here needs the coordinates.
# readr rather than data.table::fread (the idiom in diagnose_maxseg.R) so this script runs on
# a machine without data.table — a diagnostic that cannot start on a broken box is no use.
seg_sample_counts <- function(path) {
  if (!file.exists(path)) return(NULL)
  d <- tryCatch(read_tsv(path, col_select = 1, show_col_types = FALSE,
                         progress = FALSE, name_repair = "minimal"),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  as_tibble(table(d[[1]])) |> setNames(c("Sample", "n_seg"))
}

# GISTIC writes these two early, right after parsing the seg — their presence separates
# "never ingested the input" from "ingested it and crashed in the analysis".
EARLY <- c("sample_seg_counts.txt", "gistic_inputs.mat")

classify <- function(run_dir) {
  if (!dir.exists(run_dir)) return(list(state = "NO DIR", n_files = 0L))
  f <- list.files(run_dir)
  if (length(f) == 0) return(list(state = "EMPTY", n_files = 0L))
  if (any(grepl("^all_lesions\\.conf_.*\\.txt$", f)))
    return(list(state = "complete", n_files = length(f)))
  if (any(f %in% EARLY)) return(list(state = "PARTIAL", n_files = length(f)))
  list(state = "PARTIAL", n_files = length(f))
}

cat("GISTIC run diagnosis — repo ", here(), "\n",
    "  seg inputs : ", segdir, "\n",
    "  outputs    : ", outdir, "\n",
    "  -maxseg    : ", format(MAXSEG, scientific = FALSE),
    "  (export MAXSEG=... to test another cap)\n\n", sep = "")

rows <- lapply(seq_len(nrow(runs)), function(i) {
  run <- runs$run[i]
  per <- seg_sample_counts(runs$seg[i])
  cls <- classify(file.path(outdir, run))
  tibble(
    run        = run,
    seg        = if (file.exists(runs$seg[i])) "present" else "MISSING",
    n_samples  = if (is.null(per)) NA_integer_ else nrow(per),
    # THE number. GISTIC silently drops any sample over the cap, so what matters is not how
    # many samples the group has but how many it still has once the cap is applied.
    n_kept     = if (is.null(per)) NA_integer_ else sum(per$n_seg <= MAXSEG),
    max_seg    = if (is.null(per)) NA_integer_ else max(per$n_seg),
    output     = cls$state,
    n_files    = cls$n_files
  )
}) |> bind_rows()

print(as.data.frame(rows), row.names = FALSE)

cat("\n-- reading -------------------------------------------------------------\n")
bad <- rows |> filter(output != "complete")
if (nrow(bad) == 0) {
  cat("every expected run is complete.\n")
} else {
  for (i in seq_len(nrow(bad))) {
    r <- bad[i, ]
    cat("\n", r$run, ": ", r$output, "\n", sep = "")
    if (identical(r$seg, "MISSING")) {
      cat("   the input .seg was never written.\n",
          "   fix: MODE=groups Rscript code/prep_gistic_group_segs.R\n", sep = "")
    } else if (identical(r$output, "NO DIR")) {
      cat("   run_one() skipped it (no .seg at the time it ran), though one exists now.\n",
          "   fix: re-run the sweep — MODE=groups bash code/run_gistic.sh\n", sep = "")
    } else {
      # Directory present => the seg was found => re-running unchanged repeats the failure.
      cat("   the .seg was found and GISTIC itself failed; re-running as-is will fail again.\n")
      if (!is.na(r$n_kept))
        cat("   samples given: ", r$n_samples, " | surviving -maxseg ",
            format(MAXSEG, scientific = FALSE), ": ", r$n_kept,
            " | most-segmented sample: ", r$max_seg, "\n", sep = "")
      if (!is.na(r$n_kept) && r$n_kept < 5)
        cat("   ^ too few samples for GISTIC to fit a background model. Raise -maxseg so\n",
            "     fewer samples are dropped, or accept this group cannot be run alone.\n", sep = "")
    }
  }
  cat("\nRe-run only what failed by pointing OUTDIR at one run, or re-run the sweep:\n",
      "  MODE=groups bash code/run_gistic.sh\n", sep = "")
}

quit(status = if (nrow(bad) == 0) 0L else 1L)
