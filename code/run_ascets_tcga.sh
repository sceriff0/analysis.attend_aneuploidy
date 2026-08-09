#!/usr/bin/env bash
#SBATCH --job-name=ascets_tcga
#SBATCH --output=ascets_tcga_%j.log
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
# =============================================================================
# run_ascets_tcga.sh — ASCETS on the TCGA UCEC 2013 segments (report 16)
#
# WHY THIS EXISTS. ATTEND's aneuploidy score is ASCETS: the FRACTION of evaluable
# chromosome arms called AMP or DEL, i.e.
#     aneuploidy_score = #{call in AMP,DEL} / #{call != LOWCOV}          -> range 0-1
# cBioPortal's ANEUPLOIDY_SCORE is a different statistic: the Taylor et al. (Cancer Cell
# 2018) COUNT of altered arms out of 39 autosomal arms -> range 0-39, integer-valued.
# Same construct, different denominator — so the two are NOT on a common scale, and
# dividing by 39 is only approximate because ASCETS's denominator is the number of
# *evaluable* arms (LOWCOV arms drop out), which varies per sample.
#
# Running ASCETS on TCGA's own segments removes that mismatch, plus the calling-threshold
# and arm-set mismatches, leaving only the differences that are actually of interest
# (metastatic vs primary) and one that cannot be removed here (ATTEND is tumour-only
# DRAGEN WES, TCGA is purity/ploidy-corrected SNP6). Report 16 then plots ATTEND and
# TCGA on ONE scale, and report 12's standing caveat — "for an absolute-level claim,
# recompute TCGA aneuploidy with the same tool" — is discharged.
#
# INPUT   data/tcga_2013/data_cna_hg19.seg   (written by code/fetch_tcga_ucec_2013.R)
#         Already exactly ASCETS's positional format:
#           Sample  Chromosome  Start  End  Num_Probes  Segment_Mean
#         ASCETS reads these BY POSITION (ascets_resources.R renames to
#         sample/chrom/segment_start/segment_end/num_mark/log2ratio), so the header
#         names do not matter but the ORDER and the count (exactly 6) do.
#
# OUTPUT  data/tcga_2013/ascets/tcga_ucec_2013_aneuploidy_scores.txt   <- the score
#         data/tcga_2013/ascets/tcga_ucec_2013_arm_level_calls.txt     <- AMP/DEL/NEUTRAL/LOWCOV
#         data/tcga_2013/ascets/tcga_ucec_2013_arm_weighted_average_segmeans.txt
#         (+ *_params.txt, *_noise_hist.pdf, *_segmean_hist.pdf)
#         load_tcga_2013_ascets() reads the first two; report 16 prefers them over the
#         cBioPortal score automatically once they exist.
#
# RUN
#   local :  ./code/run_ascets_tcga.sh
#   HPC   :  sbatch code/run_ascets_tcga.sh
#
# ⚠️  USE THE SAME THRESHOLDS YOU USED FOR ATTEND. That is the entire point of this
#     script — a matched denominator is worthless if the calling stringency differs.
#     ✅ ALREADY MATCHED — the defaults below reproduce ATTEND's settings exactly, so a
#     plain `bash code/run_ascets_tcga.sh` is correct. ATTEND calls ascets() directly as
#         ascets(cna = cna, cytoband = cyto, min_boc = 0.5, name = sampleid,
#                alteration_threshold = 0.7)
#     which leaves TWO arguments at their defaults, and those defaults are what matter:
#       threshold = 0.2      -> amp +0.2 / del -0.2   (this script's -t 0.2)
#       noise     = data.frame()  (empty) -> ascets() takes the `else` branch and uses the
#                                 FIXED threshold; NO LCR noise modelling is performed.
#     So there is no -e to match, and the two cohorts share identical calling stringency.
#     Override MIN_BOC/ALT_FRAC/THRESHOLD only if the ATTEND run itself changes.
#
#     POOLED vs PER-SAMPLE. ATTEND runs ascets() once per sample (name = sampleid); this
#     script runs it once over all samples. That is EQUIVALENT here: the only cohort-level
#     computation in ASCETS is determine_noise_threshold(), which runs solely when a noise
#     data frame is supplied. With a fixed threshold every arm call depends on that sample's
#     own segments alone, so pooled and per-sample give identical scores. (Were -e ever
#     added, this equivalence would break and the runs would have to match in grouping too.)
#
# ENV VARS (all optional)
#   ASCETS_DIR  checkout of github.com/beroukhim-lab/ascets (auto-cloned if unset/absent)
#   BUILD       hg19 (default; the build of the TCGA 2013 segments) | hg38
#   THRESHOLD   manual log2 amp/del cutoff  (-t). Default: ASCETS's own 0.2
#   MIN_BOC     minimum arm breadth of coverage (-m). Default 0.5
#   ALT_FRAC    fraction of arm that must be altered (-a). Default 0.7
#   OUT_PREFIX  output file prefix. Default tcga_ucec_2013
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${BUILD:-hg19}"
MIN_BOC="${MIN_BOC:-0.5}"
ALT_FRAC="${ALT_FRAC:-0.7}"
THRESHOLD="${THRESHOLD:-0.2}"
OUT_PREFIX="${OUT_PREFIX:-tcga_ucec_2013}"
ASCETS_DIR="${ASCETS_DIR:-$ROOT/tools/ascets}"

SEG="$ROOT/data/tcga_2013/data_cna_hg19.seg"
OUTDIR="$ROOT/data/tcga_2013/ascets"

# --- 1. input present? ------------------------------------------------------
if [[ ! -s "$SEG" ]]; then
  echo "ERROR: $SEG not found or empty." >&2
  echo "       Run:  Rscript code/fetch_tcga_ucec_2013.R" >&2
  exit 1
fi

# --- 2. ASCETS available? ---------------------------------------------------
if [[ ! -f "$ASCETS_DIR/run_ascets.R" ]]; then
  echo "ASCETS not found at $ASCETS_DIR — cloning..."
  mkdir -p "$(dirname "$ASCETS_DIR")"
  git clone --depth 1 https://github.com/beroukhim-lab/ascets.git "$ASCETS_DIR"
fi
COORDS="$ASCETS_DIR/genomic_arm_coordinates_${BUILD}.txt"
[[ -f "$COORDS" ]] || { echo "ERROR: no arm coordinates for build '$BUILD' at $COORDS" >&2; exit 1; }

# --- 2b. slim the clone's dependencies (ASCETS_SLIM=0 to disable) ------------
# run_ascets.R and ascets_resources.R open with library(tidyverse), but ASCETS only ever
# calls dplyr and tidyr verbs (mutate/filter/group_by/summarize/left_join/arrange/distinct
# and gather/spread/separate/drop_na/replace_na). It never touches modelr, broom, ggplot2,
# forcats, readr, purrr or stringr — its plots are base graphics (pdf/hist/abline).
#
# That matters on a cluster: pulling the tidyverse META-package drags in ~8 packages ASCETS
# does not use, and ANY of them failing to build takes the whole run down. modelr failing
# to compile, or ggplot2 -> scales -> RColorBrewer being absent, both kill a run that needs
# neither. Narrowing the two library() calls removes that entire failure surface.
#
# The clone lives in tools/ (gitignored, re-clonable), so patching it is contained and
# reversible — a .orig backup is kept, and the edit is idempotent.
if [[ "${ASCETS_SLIM:-1}" == "1" ]]; then
  for f in "$ASCETS_DIR/run_ascets.R" "$ASCETS_DIR/ascets_resources.R"; do
    [[ -f "$f" ]] || continue
    if grep -q "library(tidyverse)" "$f"; then
      [[ -f "$f.orig" ]] || cp "$f" "$f.orig"
      # BRACES ARE REQUIRED: the call site is suppressMessages(suppressWarnings(library(
      # tidyverse))), so a bare `library(dplyr); library(tidyr)` would put a `;` inside a
      # function call — a syntax error. `{...}` keeps it one expression.
      sed 's/library(tidyverse)/{library(dplyr); library(tidyr)}/' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
      echo "slimmed $(basename "$f"): library(tidyverse) -> {library(dplyr); library(tidyr)}"
    fi
  done
fi

# --- 3. preflight: format + chromosome-naming compatibility ------------------
# ASCETS renames the seg columns BY POSITION, so a wrong column count silently
# mislabels everything (e.g. End read as num_mark). Fail loudly instead.
NCOL=$(awk -F'\t' 'NR==2{print NF; exit}' "$SEG")
[[ "$NCOL" == "6" ]] || { echo "ERROR: $SEG has $NCOL columns, ASCETS needs exactly 6 in the order sample/chrom/start/end/num_mark/log2ratio" >&2; exit 1; }

# The seg and the arm-coordinate table must agree on chromosome naming ("1" vs "chr1"),
# otherwise every arm silently becomes LOWCOV and every aneuploidy score comes out NaN/0.
SEG_CHR=$(awk -F'\t' 'NR==2{print $2; exit}' "$SEG")
CRD_CHR=$(awk -F'\t' 'NR==2{print $1; exit}' "$COORDS")
seg_pfx=$([[ "$SEG_CHR" == chr* ]] && echo yes || echo no)
crd_pfx=$([[ "$CRD_CHR" == chr* ]] && echo yes || echo no)
if [[ "$seg_pfx" != "$crd_pfx" ]]; then
  echo "Chromosome naming differs (seg='$SEG_CHR', coords='$CRD_CHR') — normalising the seg into a temp copy."
  SEG_FIXED="$OUTDIR/.seg_chrnorm.seg"
  mkdir -p "$OUTDIR"
  if [[ "$crd_pfx" == "yes" ]]; then
    awk -F'\t' -v OFS='\t' 'NR==1{print;next}{if($2 !~ /^chr/) $2="chr"$2; print}' "$SEG" > "$SEG_FIXED"
  else
    awk -F'\t' -v OFS='\t' 'NR==1{print;next}{sub(/^chr/,"",$2); print}' "$SEG" > "$SEG_FIXED"
  fi
  SEG="$SEG_FIXED"
fi

# --- 3b. preflight: ASCETS's own R dependencies ------------------------------
# run_ascets.R does library(tidyverse) + library(data.table) before parsing any argument,
# so a missing dependency surfaces as an R stack trace AFTER this script has printed its
# banner — which reads like the script failed. Check first and name the actual culprit.
# The usual one is RColorBrewer: `scales` IMPORTS it and ggplot2 imports scales, so an
# absent 55 KB, zero-dependency package takes the whole tidyverse down with it.
Rscript -e '
need <- if (Sys.getenv("ASCETS_SLIM", "1") == "1") c("dplyr", "tidyr", "data.table") else
        c("tidyverse", "data.table")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  cat("\nERROR: ASCETS needs these R packages and they failed to load:\n  ",
      paste(miss, collapse = ", "), "\n", sep = "")
  # A top-level package can fail because a TRANSITIVE dependency is absent; name it.
  deep <- c("RColorBrewer", "scales", "ggplot2", "dplyr", "tidyr", "readr", "stringr", "purrr")
  gone <- deep[!vapply(deep, function(p) length(find.package(p, quiet = TRUE)) > 0, logical(1))]
  if (length(gone))
    cat("Actually missing from the library path:\n  ", paste(gone, collapse = ", "), "\n",
        "Install with:\n  install.packages(c(\"", paste(gone, collapse = "\",\""),
        "\"), repos = \"https://cloud.r-project.org\")\n", sep = "")
  cat("Library paths searched:\n  ", paste(.libPaths(), collapse = "\n  "), "\n", sep = "")
  quit(status = 1)
}' || exit 1

# --- 4. run -----------------------------------------------------------------
# run_ascets.R resolves ascets_resources.R from its OWN path, but writes its outputs to
# the CURRENT WORKING DIRECTORY (write_outputs_to_file() defaults location="./"), so we
# cd into the output folder and pass absolute paths in.
mkdir -p "$OUTDIR"
echo "=== ASCETS on TCGA UCEC 2013 ==="
echo "  seg        : $SEG"
echo "  coords     : $COORDS  (build $BUILD)"
echo "  thresholds : -t $THRESHOLD  -m $MIN_BOC  -a $ALT_FRAC"
echo "  out        : $OUTDIR/${OUT_PREFIX}_*"
cd "$OUTDIR"
Rscript "$ASCETS_DIR/run_ascets.R" \
  -i "$SEG" \
  -c "$COORDS" \
  -m "$MIN_BOC" \
  -a "$ALT_FRAC" \
  -t "$THRESHOLD" \
  -o "$OUT_PREFIX"

rm -f "$OUTDIR/.seg_chrnorm.seg"

# --- 5. summarise -----------------------------------------------------------
SCORES="$OUTDIR/${OUT_PREFIX}_aneuploidy_scores.txt"
if [[ -s "$SCORES" ]]; then
  echo
  awk -F'\t' 'NR>1{n++; s+=$2; a[n]=$2} END{
    if(n==0){print "  no scores written!"; exit}
    asort(a); m=(n%2)?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2
    printf "  %d samples scored; mean %.4f, median %.4f (range 0-1)\n", n, s/n, m
  }' "$SCORES" 2>/dev/null || echo "  $(( $(wc -l < "$SCORES") - 1 )) samples scored -> $SCORES"
  echo
  echo "Next: knit analysis/16_tcga2013_fig1a.Rmd — it prefers this score automatically."
else
  echo "WARNING: expected $SCORES but it is missing/empty." >&2
  exit 1
fi
