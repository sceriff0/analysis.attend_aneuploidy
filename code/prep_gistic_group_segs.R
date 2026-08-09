#!/usr/bin/env Rscript
# =============================================================================
# prep_gistic_group_segs.R — write the per-GROUP and leave-one-out .seg inputs
# that code/run_gistic.sh consumes in MODE=groups / MODE=loo.
#
# WHY THIS EXISTS: write_gistic_seg() alone emits ONLY the pooled seg
# (output/gistic_input/attend_all_segments.seg). With no attend_<group>.seg /
# attend_mmrd_high_drop_*.seg present, run_one() in run_gistic.sh hits its
# `skip: <file> missing` guard and MODE=groups/loo silently produce nothing —
# which is why only data/gistic/all/ ever appeared.
#
# It also writes the POOLED seg (attend_all_segments.seg) itself, so you do NOT need
# the `Rscript -e 'source(...); write_gistic_seg()'` one-liner (which fails with
# "object 'attend_cnv' not found" — that command never sources attend_classes.R).
#
# RUN (from the repo root, on the cluster where data/seg/ lives). MODE mirrors
# run_gistic.sh's MODE, so prep and GISTIC take the SAME value:
#   Rscript code/prep_gistic_group_segs.R              # pooled + groups + LOO (default: both)
#   MODE=pooled Rscript code/prep_gistic_group_segs.R  # pooled seg only
#   MODE=groups Rscript code/prep_gistic_group_segs.R  # pooled + 4 group segs
#   MODE=loo    Rscript code/prep_gistic_group_segs.R  # MMRd-high LOO segs only
#
# THEN feed GISTIC (each mode into its OWN namespaced dir — no overwrite):
#   MODE=groups bash code/run_gistic.sh                # pooled all/ + the 4 group runs
#   MODE=loo    bash code/run_gistic.sh                # the LOO runs
#
# Prereqs: report 1 has written the master; data/seg/*.seg present;
#          load_gianlu_clinical_data() resolves (barcode->pid crosswalk).
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(data.table); library(fs)
})
source(here("code", "attend_classes.R"))
source(here("code", "attend_scna.R"))
source(here("code", "attend_harmonise.R"))   # make_id_cfg(), build_barcode_pid()
source(here("code", "attend_io.R"))          # read_intermediate()
source(here("code", "load_wes_results.R"))   # write_gistic_seg()
source(here("code", "load_clinical.R"))      # load_gianlu_clinical_data()

mode <- Sys.getenv("MODE", "both")
stopifnot("MODE must be pooled | groups | loo | both" =
            mode %in% c("pooled", "groups", "loo", "both"))
need_pooled <- mode %in% c("pooled", "groups", "both")   # groups' all/ run needs it
need_groups <- mode %in% c("groups", "both")
need_loo    <- mode %in% c("loo",    "both")

# --- pooled seg (attend_all_segments.seg) — needs NO master/crosswalk ---------
if (need_pooled) {
  cat(">> pooled seg\n")
  if (is.null(write_gistic_seg()))
    stop("write_gistic_seg() wrote nothing — is data/seg/ populated? (see attend_cnv$seg$dir)")
}

# groups/LOO need the grouping; pooled-only skips this whole block.
if (need_groups || need_loo) {
  # --- master + scna_group (same composition report 5 uses) -----------------
  master <- tryCatch(read_intermediate("attend_master_joined"), error = function(e) NULL)
  if (is.null(master) || !nrow(master))
    stop("master table absent — knit report 1 first (output/clean_data/attend_master_joined).")
  master <- add_scna_group(add_molecular_classes(master))

  # --- barcode -> pid crosswalk (same derivation every other report/15 use) --
  cw <- tryCatch(build_barcode_pid(load_gianlu_clinical_data(), make_id_cfg()),
                 error = function(e) NULL)
  if (is.null(cw) || !nrow(cw))
    stop("barcode->pid crosswalk empty — check load_gianlu_clinical_data() / make_id_cfg().")

  # --- sanity: group sizes among the actual .seg barcodes, BEFORE GISTIC runs -
  # A zero/tiny group here is far cheaper to catch now than after a multi-hour run.
  bc  <- seg_barcodes()
  grp <- barcode_scna_group(bc, master, cw)
  cat("scna_group membership among", length(bc), "seg-file barcodes:\n")
  print(table(scna_group = factor(as.character(grp), levels = attend_scna$group_levels),
              useNA = "ifany"))
  cat("\n")

  if (need_groups) write_group_segs(master, cw)
  if (need_loo)    write_loo_segs(master, cw)
}

cat("\nNext (each mode writes its OWN data/gistic/<...> dir — no overwrite):\n",
    "  MODE=groups bash code/run_gistic.sh\n",
    "  MODE=loo    bash code/run_gistic.sh\n", sep = "")
