#!/usr/bin/env bash
#SBATCH --job-name=gistic_attend
#SBATCH --output=gistic_%j.log
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
# =============================================================================
# run_gistic.sh — GISTIC2.0 on the ATTEND DRAGEN copy-number segments
#
# GISTIC clusters copy-number profiles for report 15 (recurrent SCNA by aneuploidy/MMR).
# Earlier: report 09 used GISTIC output for TCGA integration. It runs OUTSIDE R (a
# compiled MATLAB binary) on per-group seg files R writes via write_gistic_seg().
# Outputs land in data/gistic/<group>/, where report 15 and find_gistic_files()
# (report 06 oncoplot) pick them up.
#
# Runs BOTH ways:
#   local :  ./code/run_gistic.sh
#   HPC   :  sbatch code/run_gistic.sh     (the #SBATCH lines above are used)
#
# FULL PIPELINE (on the cluster where the .seg live):
#   1. sync/keep the DRAGEN .seg in data/seg/  (dragen_standard)
#   2. R:   Rscript code/prep_gistic_group_segs.R   (MODE=pooled|groups|loo|both)
#           writes the pooled seg AND the per-group / LOO segs -> output/gistic_input/
#           (do NOT use `Rscript -e 'source("code/load_wes_results.R"); write_gistic_seg()'`
#            — it never sources attend_classes.R, so attend_cnv is undefined and it errors.)
#   3. bash code/run_gistic.sh  or  MODE=groups bash code/run_gistic.sh  or  MODE=loo bash code/run_gistic.sh
#      (edit GISTIC_SIF / REFGENE / MODULE below)
#   4. R:   knit analysis/15_recurrent_scna_by_aneuploidy_mmr.Rmd
#
# The DRAGEN .seg is ALREADY GISTIC's format (Sample/Chr/Start/End/Num_Markers/Seg.CN),
# and seg.mean is log2 — no transform needed.
#
# GETTING GISTIC ON HPC (pick the one your cluster supports):
#   A. Singularity/Apptainer (most portable):
#        singularity pull gistic2.sif docker://genepattern/docker-gistic
#        GISTIC_SIF=gistic2.sif        # this script then wraps it automatically
#   B. environment module:   module load gistic/2.0.23      # sets `gistic2` on PATH
#   C. conda:                conda activate gistic2         # bioconda gistic2 (bundles MCR)
#
# ⚠️ hg38 REFERENCE (critical). The stock binary ships hg16-19 only; ATTEND is hg38.
#   Supply hg38.UCSC.add_mir.160920.refgene.mat (GISTIC2 GitHub / bzhanglab/GISTIC2_example).
#   Running hg19 on hg38 segments misplaces EVERY peak.
# =============================================================================
source ~/.bashrc
cd "${SLURM_SUBMIT_DIR:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"

# --- edit for your environment ----------------------------------------------
MODULE="${MODULE:-}"                    # e.g. "gistic/2.0.23"; leave empty if not using modules
GISTIC_SIF="${GISTIC_SIF:-}"            # path to a Singularity image; if set, we wrap it
GISTIC_BIN="${GISTIC_BIN:-gistic2}"     # used when GISTIC_SIF is empty (module/conda/PATH)
REFGENE="${REFGENE:-refgenefiles/hg38.UCSC.add_mir.160920.refgene.mat}"
SEG="${SEG:-output/gistic_input/attend_all_segments.seg}"
OUTDIR="${OUTDIR:-data/gistic}"
# GISTIC EXCLUDES any sample whose segment count exceeds -maxseg (silently — it just
# vanishes from all_thresholded.by_genes.txt). DRAGEN WES .seg are finely segmented
# (often >>2000 segments/sample). ATTEND uses 46000 as a deliberate hypersegmentation
# filter for report 15's small MMRd-high group — excluded samples are audited in
# report 15 via maxseg_audit() (see: Rscript code/diagnose_maxseg.R).
MAXSEG="${MAXSEG:-46000}"
# ----------------------------------------------------------------------------

[ -n "$MODULE" ] && { module load "$MODULE" 2>/dev/null || echo "note: 'module load $MODULE' failed — continuing"; }
if [ -n "$GISTIC_SIF" ]; then
  # $GISTIC_BIN is the command name/path INSIDE the container. Some images (e.g.
  # genepattern/docker-gistic) don't put `gistic2` on PATH — set GISTIC_BIN to the
  # full in-container path (e.g. /usr/local/bin/gp_gistic2_from_seg) in that case.
  RUN=(singularity exec --bind "$PWD:$PWD" --pwd "$PWD" "$GISTIC_SIF" "$GISTIC_BIN")
else
  RUN=("$GISTIC_BIN")
fi

# MODE=pooled  : the pooled run only (default)
# MODE=groups  : pooled + the four scna_group runs
# MODE=loo     : the 9 leave-one-out runs for MMRd-high (peak stability, spec §4.2)
MODE="${MODE:-pooled}"
SEGDIR="${SEGDIR:-output/gistic_input}"

# $SEG is only needed by pooled/groups (the pooled run); loo uses $SEGDIR instead.
case "$MODE" in
  pooled|groups) [ -f "$SEG" ] || { echo "ERROR: $SEG missing — run write_gistic_seg() in R first."; exit 1; } ;;
esac
[ -f "$REFGENE" ] || echo "WARNING: refgene '$REFGENE' not found — set REFGENE to your hg38 .mat"
mkdir -p "$OUTDIR"

# Pre-registered parameters (see docs/superpowers/specs/2026-07-20-recurrent-scna-by-aneuploidy-mmr-design.md §4.3).
# Stricter -ta/-td 0.3 and -conf 0.99 give fewer, narrower, higher-amplitude peaks —
# the correct FDR posture for report 15's small MMRd-high group.
run_one () {
  local seg_file="$1" out_dir="$2"
  [ -f "$seg_file" ] || { echo "skip: $seg_file missing"; return 0; }
  mkdir -p "$out_dir"
  "${RUN[@]}" \
    -b "$out_dir" \
    -seg "$seg_file" \
    -refgene "$REFGENE" \
    -genegistic 1 \
    -broad 1 -brlen 0.7 \
    -conf 0.99 \
    -armpeel 1 \
    -savegene 1 \
    -gcm extreme \
    -smallmem 0 \
    -rx 0 \
    -ta 0.3 -td 0.3 \
    -cap 1.5 \
    -v 30 \
    -js 4 \
    -maxseg "$MAXSEG"
  echo "  -> $out_dir"
}

case "$MODE" in
  pooled)
    run_one "$SEG" "$OUTDIR/all"
    ;;
  groups)
    run_one "$SEG" "$OUTDIR/all"
    # || echo WARN keeps the sweep going if one run fails (set -e is suppressed in a
    # function called in a || list) — one failing group must not abort the other three.
    for g in mmrp_high mmrp_low mmrd_low mmrd_high; do
      echo "GISTIC group: $g"
      run_one "$SEGDIR/attend_${g}.seg" "$OUTDIR/$g" || echo "WARN: GISTIC failed for group $g — continuing"
    done
    ;;
  loo)
    # Leave-one-out over MMRd-high. R writes one attend_mmrd_high_drop_<id>.seg per
    # sample via write_gistic_seg(ids = setdiff(mmrd_high_ids, id)).
    # || echo WARN keeps the sweep going if one run fails (set -e is suppressed in a
    # function called in a || list) — LOO requires all 9 runs to ATTEMPT.
    for f in "$SEGDIR"/attend_mmrd_high_drop_*.seg; do
      [ -e "$f" ] || { echo "no LOO seg files in $SEGDIR"; break; }
      b=$(basename "$f" .seg); b=${b#attend_}
      echo "GISTIC LOO: $b"
      run_one "$f" "$OUTDIR/loo/$b" || echo "WARN: GISTIC failed for LOO $b — continuing"
    done
    ;;
  *)
    echo "ERROR: unknown MODE='$MODE' (expected pooled|groups|loo)"; exit 1
    ;;
esac

echo "GISTIC2 done (MODE=$MODE) -> $OUTDIR"
echo "  clustering input : $OUTDIR/all/all_thresholded.by_genes.txt"
echo "  peak calls       : <group>/all_lesions.conf_99.txt  (report 15 -> load_gistic_lesions_at)"
echo "Next: knit analysis/15_recurrent_scna_by_aneuploidy_mmr.Rmd"
