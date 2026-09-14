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
#   Supply hg38.UCSC.add_miR.160920.refgene.mat (GISTIC2 GitHub / bzhanglab/GISTIC2_example).
#   Running hg19 on hg38 segments misplaces EVERY peak.
# =============================================================================
# --- locate the repo root, and PROVE it ---------------------------------------------------
# The old line was `cd "${SLURM_SUBMIT_DIR:-$(git rev-parse --show-toplevel || pwd)}"`, which
# fails three ways on a cluster: sbatch sets SLURM_SUBMIT_DIR to wherever you INVOKED sbatch
# (not the repo), `git rev-parse` needs the cwd to already be inside the repo, and the final
# `pwd` fallback silently accepts any directory at all — after which every relative path
# below resolves somewhere wrong and GISTIC writes its output into the void.
#
# Candidates are tried in order and each is VALIDATED against a sentinel (code/run_gistic.sh
# + analysis/ must both exist) before being accepted, so a wrong guess is rejected rather
# than used. If none validates the script says which ones it tried and stops.
#
# NOTE ON $0 UNDER SLURM: sbatch copies the batch script into its spool directory, so
# BASH_SOURCE can point at /var/spool/.../slurm_script. That is why the script-location
# candidate is validated like every other one rather than trusted.
_valid_root () { [ -f "$1/code/run_gistic.sh" ] && [ -d "$1/analysis" ]; }

_self="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" 2>/dev/null && pwd)"
_cands=(
  "${ATTEND_ROOT:-}"                                   # explicit override, wins outright
  "${_self%/code}"                                     # the script's own repo, when not spooled
  "${SLURM_SUBMIT_DIR:-}"                              # where sbatch was invoked
  "$(git rev-parse --show-toplevel 2>/dev/null || true)"
  "$PWD"
  "$HOME/workflowR/attend_aneuploidy"                  # the cluster checkout
)
ROOT=""
for _c in "${_cands[@]}"; do
  [ -n "$_c" ] || continue
  if _valid_root "$_c"; then ROOT="$(cd "$_c" && pwd)"; break; fi
done
if [ -z "$ROOT" ]; then
  echo "ERROR: could not find the attend_aneuploidy repo root."
  echo "       A valid root contains BOTH code/run_gistic.sh and analysis/."
  echo "       Tried, in order:"
  for _c in "${_cands[@]}"; do [ -n "$_c" ] && echo "         - $_c"; done
  echo "       Fix: ATTEND_ROOT=~/workflowR/attend_aneuploidy sbatch code/run_gistic.sh"
  exit 1
fi
cd "$ROOT" || exit 1
echo "repo root: $ROOT"

# ~/.bashrc is sourced ONLY for `module`, and only when it is safe to. Most .bashrc files
# begin with "if not running interactively, return", so under sbatch this is frequently a
# no-op — which is why `module load` can work on the login node and silently not work in the
# job. If `module` is still undefined afterwards the MODULE branch below says so rather than
# failing with "command not found" halfway through.
if [ -n "${MODULE:-}" ] && ! command -v module >/dev/null 2>&1; then
  # shellcheck disable=SC1090
  [ -f "$HOME/.bashrc" ] && source "$HOME/.bashrc" 2>/dev/null || true
  command -v module >/dev/null 2>&1 || \
    echo "note: MODULE='$MODULE' requested but \`module\` is not available in this shell." \
         "On many clusters ~/.bashrc returns early when non-interactive — source the" \
         "modules init directly (e.g. /etc/profile.d/modules.sh) or use GISTIC_SIF instead."
fi

# --- edit for your environment ----------------------------------------------
MODULE="${MODULE:-}"                    # e.g. "gistic/2.0.23"; leave empty if not using modules
GISTIC_SIF="${GISTIC_SIF:-}"            # path to a Singularity image; if set, we wrap it
GISTIC_BIN="${GISTIC_BIN:-gistic2}"     # used when GISTIC_SIF is empty (module/conda/PATH)
REFGENE="${REFGENE:-refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat}"
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

# Every path is made absolute against $ROOT. A relative REFGENE/SEG/OUTDIR passed in the
# environment then means the same thing wherever the job was submitted from, and — the part
# that actually matters — the container gets real paths it can bind.
_abs () { case "$1" in /*) printf '%s' "$1" ;; *) printf '%s/%s' "$ROOT" "$1" ;; esac; }
REFGENE="$(_abs "$REFGENE")"; SEG="$(_abs "$SEG")"; OUTDIR="$(_abs "$OUTDIR")"

# GISTIC_SIF is TRIMMED, then VALIDATED, then auto-detected — in that order.
#
# Trimmed because an exported variable picks up whitespace easily and `[ -n " " ]` is true,
# so a blank-looking value silently takes the container branch. Validated because a stale
# export from an earlier shell is the likeliest way this goes wrong: a GISTIC_SIF pointing at
# a DIRECTORY (the repo root, say) passed `[ -n ]`, survived _abs() unchanged because it was
# already absolute, and died as "does not exist" naming a path that plainly does exist —
# which reads as a bug in the script rather than a stale variable. Auto-detected because the
# image normally sits in the repo root as gistic2.sif, exactly as the header instructs, so
# the common case should need no environment at all.
GISTIC_SIF="$(printf '%s' "$GISTIC_SIF" | tr -d '[:space:]')"
if [ -n "$GISTIC_SIF" ]; then
  GISTIC_SIF="$(_abs "$GISTIC_SIF")"
  if [ -d "$GISTIC_SIF" ]; then
    echo "note: GISTIC_SIF='$GISTIC_SIF' is a DIRECTORY, not a .sif — ignoring it."
    echo "      (a stale 'export GISTIC_SIF=...' in this shell? check with: echo \"[\$GISTIC_SIF]\")"
    GISTIC_SIF=""
  elif [ ! -f "$GISTIC_SIF" ]; then
    echo "note: GISTIC_SIF='$GISTIC_SIF' does not exist — ignoring it and looking in the repo."
    GISTIC_SIF=""
  fi
fi
if [ -z "$GISTIC_SIF" ]; then
  for _c in "$ROOT/gistic2.sif" "$ROOT"/*.sif; do
    [ -f "$_c" ] && { GISTIC_SIF="$_c"; echo "found container: $GISTIC_SIF"; break; }
  done
fi

if [ -n "$GISTIC_SIF" ]; then
  # ⚠️ BIND AS LITTLE AS POSSIBLE. genepattern/docker-gistic is a MINIMAL image — it has no
  # /tmp, no /var/tmp and no /etc/passwd, which singularity announces as "Skipping mount
  # /tmp: /tmp doesn't exist in container". An image that thin cannot have mount points
  # created inside it, so an explicit --bind (and a --pwd onto a bound path) can leave the
  # container in a state where singularity then fails to stat the binary itself:
  #
  #   FATAL: stat /usr/local/bin/gp_gistic2_from_seg: no such file or directory
  #
  # ...for a binary that `singularity exec gistic2.sif ls -l` shows is plainly there. The bare
  # invocation works because singularity ALREADY auto-binds $HOME and the current directory,
  # and the repo lives under $HOME. So bind nothing by default, and add a path only when it
  # is genuinely outside both — which is the only case the auto-binds miss.
  _needs_bind () {   # true when $1 is under neither $HOME nor $ROOT
    case "$1" in "$HOME"/*|"$HOME"|"$ROOT"/*|"$ROOT") return 1 ;; *) return 0 ;; esac
  }
  BIND_ARGS=(); _bound=()
  for _b in "$ROOT" "$(dirname "$GISTIC_SIF")" "$(dirname "$REFGENE")"; do
    [ -d "$_b" ] || continue
    _needs_bind "$_b" || continue
    case " ${_bound[*]} " in *" $_b "*) continue ;; esac
    _bound+=("$_b"); BIND_ARGS+=(--bind "$_b:$_b")
  done
  # MATLAB's MCR needs somewhere writable for its component cache, and this image has no
  # /tmp. Left unset it can fail deep inside GISTIC with an unhelpful MATLAB error, so it is
  # pointed at the repo unless the caller overrides it.
  export MCR_CACHE_ROOT="${MCR_CACHE_ROOT:-$ROOT/.mcr_cache}"
  mkdir -p "$MCR_CACHE_ROOT"
  # $GISTIC_BIN is the command name/path INSIDE the container. This image does NOT put
  # `gistic2` on PATH — it ships /usr/local/bin/gp_gistic2_from_seg. List what is there with:
  #   singularity exec gistic2.sif bash -lc 'ls /usr/local/bin /opt/GISTIC 2>/dev/null'
  RUN=(singularity exec "${BIND_ARGS[@]}" "$GISTIC_SIF" "$GISTIC_BIN")
  echo "gistic: $GISTIC_SIF -> $GISTIC_BIN"
  echo "  binds: ${_bound[*]:-none needed (\$HOME and cwd are auto-bound)}"
  echo "  MCR_CACHE_ROOT: $MCR_CACHE_ROOT"
  # Prove the binary is reachable BEFORE launching five multi-hour runs against it, and find
  # it if the configured name is wrong. The default is `gistic2`, which this image does NOT
  # ship — genepattern/docker-gistic puts /usr/local/bin/gp_gistic2_from_seg instead — so
  # without the probe the zero-config path fails for a reason the FATAL does not name.
  _bin_ok () { singularity exec "${BIND_ARGS[@]}" "$GISTIC_SIF" test -x "$1" 2>/dev/null; }
  if ! _bin_ok "$GISTIC_BIN"; then
    _found=""
    for _b in /usr/local/bin/gp_gistic2_from_seg /opt/GISTIC/gp_gistic2_from_seg \
              /opt/GISTIC/gistic2 gp_gistic2_from_seg gistic2; do
      _bin_ok "$_b" && { _found="$_b"; break; }
    done
    if [ -n "$_found" ]; then
      echo "note: GISTIC_BIN='$GISTIC_BIN' is not executable in the image; using '$_found'."
      GISTIC_BIN="$_found"
    else
      echo "ERROR: no GISTIC binary found inside $GISTIC_SIF."
      echo "       Tried: $GISTIC_BIN, /usr/local/bin/gp_gistic2_from_seg,"
      echo "              /opt/GISTIC/gp_gistic2_from_seg, /opt/GISTIC/gistic2, gistic2"
      echo "       List what is actually there with:"
      echo "         singularity exec $GISTIC_SIF bash -lc 'ls -l /usr/local/bin /opt/GISTIC 2>/dev/null'"
      exit 1
    fi
  fi
  RUN=(singularity exec "${BIND_ARGS[@]}" "$GISTIC_SIF" "$GISTIC_BIN")
  echo "  binary: $GISTIC_BIN (verified executable in the image)"
else
  command -v "$GISTIC_BIN" >/dev/null 2>&1 || {
    echo "ERROR: '$GISTIC_BIN' is not on PATH and GISTIC_SIF is unset."
    echo "       Set one of: GISTIC_SIF=<image>.sif | MODULE=gistic/2.0.23 | GISTIC_BIN=<path>"
    exit 1; }
  RUN=("$GISTIC_BIN")
  echo "gistic: $(command -v "$GISTIC_BIN")"
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
FAILED_RUNS=()   # names of runs that did not produce output; summarised at the end

run_one () {
  local seg_file="$1" out_dir="$2" rc=0
  [ -f "$seg_file" ] || { echo "skip: $seg_file missing"; FAILED_RUNS+=("$(basename "$out_dir") [no .seg]"); return 0; }
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
    -maxseg "$MAXSEG" || rc=$?

  # ⚠️ A RUN THAT FAILED MUST NOT REPORT SUCCESS. This used to `echo "  -> $out_dir"`
  # unconditionally, so a singularity FATAL or a GISTIC crash printed an arrow and the script
  # went on to announce "GISTIC2 done". mkdir -p had already created the folder, so the only
  # trace was an output directory holding nothing — indistinguishable on disk from a run that
  # was never attempted, and the reports downstream just found no peaks.
  #
  # Exit status alone is not enough: GISTIC can exit 0 having written nothing. So the check is
  # BOTH — non-zero status, or the absence of the file every consumer needs.
  local lesions; lesions=$(ls "$out_dir"/all_lesions.conf_*.txt 2>/dev/null | head -1)
  if [ "$rc" -ne 0 ] || [ -z "$lesions" ]; then
    echo "  !! FAILED: $out_dir  (exit $rc; all_lesions.conf_*.txt $( [ -n "$lesions" ] && echo present || echo absent ))"
    FAILED_RUNS+=("$(basename "$out_dir")")
    return 0          # keep the sweep going; the summary below is what reports the failures
  fi
  echo "  -> $out_dir  (ok)"
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

if [ "${#FAILED_RUNS[@]}" -gt 0 ]; then
  echo
  echo "=============================================================================="
  echo "GISTIC2 (MODE=$MODE): ${#FAILED_RUNS[@]} run(s) produced NO OUTPUT:"
  for r in "${FAILED_RUNS[@]}"; do echo "  - $r"; done
  echo
  echo "Those directories exist but are empty or partial, which on disk looks exactly"
  echo "like a run that was never attempted. Reports reading them will find no peaks."
  echo "Check the FATAL/error lines above; common causes:"
  echo "  - singularity cannot stat the in-container binary -> GISTIC_BIN wrong, or a"
  echo "    --bind is shadowing the image. Verify with:"
  echo "      singularity exec $GISTIC_SIF ls -l $GISTIC_BIN"
  echo "  - MATLAB MCR has nowhere to write -> export MCR_CACHE_ROOT=/path/on/scratch"
  echo "  - refgene build mismatch -> REFGENE must be the hg38 .mat"
  echo "=============================================================================="
  exit 1
fi

echo "GISTIC2 done (MODE=$MODE) -> $OUTDIR   (all runs produced output)"
echo "  clustering input : $OUTDIR/all/all_thresholded.by_genes.txt"
echo "  peak calls       : <group>/all_lesions.conf_99.txt  (-> load_gistic_lesions_at)"
echo "Next: knit analysis/10-recurrent-cna.Rmd  (and 07-molecular-classification.Rmd)"
