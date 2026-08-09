#!/usr/bin/env bash
#SBATCH --job-name=gistic_tcga
#SBATCH --output=gistic_tcga_%j.log
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
# =============================================================================
# run_gistic_tcga.sh — GISTIC2.0 on the TCGA UCEC 2013 segments (report 16)
#
# Sibling of code/run_gistic.sh, which runs GISTIC on the ATTEND DRAGEN segments. This
# one runs the SAME pre-registered parameters on TCGA's own SNP6 segments so report 16
# can cluster the TCGA cohort on the *faithful* TCGA feature space — Kandoth et al.
# Suppl. Methods S2: "thresholded relative copy number data in significantly reoccurring
# amplification or deletion regions identified by GISTIC 2.0" — instead of the naive
# genome-bin matrix. Report 16 then shows four heatmaps of one cohort:
#     1. published CNA_CLUSTER_K4        (TCGA's own answer)
#     2. genome-bin clustering           (naive; no GISTIC)
#     3. GISTIC peaks, Euclidean + Ward
#     4. GISTIC peaks, Pearson + Ward
# and scores 2-4 against 1 with an adjusted Rand index.
#
# INPUT   data/tcga_2013/data_cna_hg19.seg   (from code/fetch_tcga_ucec_2013.R)
#         Already GISTIC's format: Sample / Chromosome / Start / End / Num_Probes /
#         Segment_Mean, and Segment_Mean is log2 — no transform needed, same as DRAGEN.
# OUTPUT  data/tcga_2013/gistic/all/   — all_thresholded.by_genes.txt (the clustering
#         input), amp_genes/del_genes (the peak restriction), all_lesions, scores.gistic.
#         Report 16 reads them via load_gistic_thresholded(cnv = attend_tcga_gistic).
#
# RUN
#   local :  ./code/run_gistic_tcga.sh
#   HPC   :  sbatch code/run_gistic_tcga.sh
#
# ============================ TWO DIFFERENCES FROM run_gistic.sh =============
#
# 1. ✅ REFERENCE GENOME IS EASY HERE. run_gistic.sh carries a loud warning because
#    ATTEND is hg38 and the stock GISTIC binary ships hg16-hg19 only, so it needs a
#    hand-supplied hg38 refgene .mat. TCGA 2013 is **hg19**, which the binary ships
#    natively — so this script just picks the bundled hg19 file up. Nothing to source.
#
# 2. ✅ -maxseg IS NOT A FILTER HERE. ATTEND sets -maxseg 46000 as a deliberate
#    hypersegmentation filter: DRAGEN WES segments are fine-grained (often >>2000
#    segments/sample) and GISTIC SILENTLY DROPS any sample above the limit. TCGA's SNP6
#    segments are coarse — 143 per sample on average, 1966 at the maximum — so the
#    default 10000 below is pure headroom and excludes nobody. Report 16 checks the
#    surviving sample count against the seg file regardless, so a silent drop cannot
#    pass unnoticed.
#
# ⚠️ WHAT DOES *NOT* TRANSFER: the peak calls themselves. -ta/-td 0.3 and -conf 0.99 were
#    pre-registered for report 15's small MMRd-high group (fewer, narrower, higher-
#    amplitude peaks). They are kept here so the ATTEND and TCGA feature spaces are built
#    the same way — which is the point — but the resulting TCGA peak list is NOT the one
#    published in 2013 (different GISTIC version, parameters and reference). Report 16
#    therefore compares CLUSTERINGS, never peak lists.
#
# GETTING GISTIC (same options as run_gistic.sh):
#   A. singularity pull gistic2.sif docker://genepattern/docker-gistic ; GISTIC_SIF=gistic2.sif
#   B. module load gistic/2.0.23
#   C. conda activate gistic2
# =============================================================================
source ~/.bashrc 2>/dev/null || true
cd "${SLURM_SUBMIT_DIR:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"

# --- edit for your environment ----------------------------------------------
MODULE="${MODULE:-}"                    # e.g. "gistic/2.0.23"; empty if not using modules
GISTIC_SIF="${GISTIC_SIF:-}"            # path to a Singularity image; if set, we wrap it
GISTIC_BIN="${GISTIC_BIN:-gistic2}"     # command when GISTIC_SIF is empty
SEG="${SEG:-data/tcga_2013/data_cna_hg19.seg}"
OUTDIR="${OUTDIR:-data/tcga_2013/gistic}"
MAXSEG="${MAXSEG:-10000}"               # headroom, not a filter — see note 2 above
# hg19 refgene: set REFGENE explicitly, else we try the usual bundled names in order.
REFGENE="${REFGENE:-}"
# ----------------------------------------------------------------------------

[ -n "$MODULE" ] && { module load "$MODULE" 2>/dev/null || echo "note: 'module load $MODULE' failed — continuing"; }
if [ -n "$GISTIC_SIF" ]; then
  RUN=(singularity exec --bind "$PWD:$PWD" --pwd "$PWD" "$GISTIC_SIF" "$GISTIC_BIN")
else
  RUN=("$GISTIC_BIN")
fi

[ -f "$SEG" ] || { echo "ERROR: $SEG missing — run 'Rscript code/fetch_tcga_ucec_2013.R' first."; exit 1; }

# --- locate the hg19 refgene ------------------------------------------------
# Distributions disagree on the filename; try the common ones rather than hard-coding.
if [ -z "$REFGENE" ]; then
  for cand in \
      refgenefiles/hg19.UCSC.add_mir.140312.refgene.mat \
      refgenefiles/hg19.mat \
      /opt/GISTIC/refgenefiles/hg19.UCSC.add_mir.140312.refgene.mat \
      /opt/GISTIC/refgenefiles/hg19.mat \
      "$GISTIC_REFDIR/hg19.mat"; do
    [ -f "$cand" ] && { REFGENE="$cand"; break; }
  done
fi
if [ -z "$REFGENE" ] || [ ! -f "$REFGENE" ]; then
  echo "WARNING: no hg19 refgene .mat found. Set REFGENE=/path/to/hg19...refgene.mat"
  echo "         (inside a container it is usually /opt/GISTIC/refgenefiles/hg19.mat)"
  echo "         ⚠️  These segments are hg19 — do NOT substitute an hg38 reference."
fi

# --- preflight --------------------------------------------------------------
NCOL=$(awk -F'\t' 'NR==2{print NF; exit}' "$SEG")
[ "$NCOL" = "6" ] || { echo "ERROR: $SEG has $NCOL columns; GISTIC needs 6 (Sample/Chr/Start/End/Num_Markers/Seg.CN)"; exit 1; }
awk -F'\t' 'NR>1{c[$1]++} END{n=0;mx=0;s=0; for(k in c){n++;s+=c[k]; if(c[k]>mx)mx=c[k]}
  printf "  %d samples, %d segments (mean %.0f/sample, max %d)\n", n, s, s/n, mx}' "$SEG"
MAXOBS=$(awk -F'\t' 'NR>1{c[$1]++} END{mx=0; for(k in c) if(c[k]>mx)mx=c[k]; print mx}' "$SEG")
if [ "$MAXOBS" -gt "$MAXSEG" ]; then
  echo "ERROR: a sample has $MAXOBS segments but -maxseg is $MAXSEG — GISTIC would SILENTLY DROP it."
  echo "       Raise MAXSEG above $MAXOBS."; exit 1
fi

mkdir -p "$OUTDIR/all"
echo "=== GISTIC2 on TCGA UCEC 2013 ==="
echo "  seg     : $SEG"
echo "  refgene : ${REFGENE:-<unset>}  (hg19)"
echo "  out     : $OUTDIR/all"

# Same pre-registered parameters as code/run_gistic.sh, so the ATTEND and TCGA feature
# spaces are constructed identically (docs/superpowers/specs/...-recurrent-scna-... §4.3).
"${RUN[@]}" \
  -b "$OUTDIR/all" \
  -seg "$SEG" \
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

# --- verify the sample count survived ---------------------------------------
THR="$OUTDIR/all/all_thresholded.by_genes.txt"
if [ -s "$THR" ]; then
  NIN=$(awk -F'\t' 'NR>1{c[$1]++} END{n=0; for(k in c)n++; print n}' "$SEG")
  NOUT=$(head -1 "$THR" | awk -F'\t' '{print NF-3}')
  echo
  echo "  samples in  : $NIN"
  echo "  samples out : $NOUT"
  [ "$NIN" = "$NOUT" ] || echo "  ⚠️  $((NIN-NOUT)) sample(s) dropped by GISTIC — check -maxseg and the log above."
  echo
  echo "GISTIC2 done -> $OUTDIR/all"
  echo "  clustering input : $THR"
  echo "Next: knit analysis/16_tcga2013_fig1a.Rmd (it picks these up automatically)."
else
  echo "ERROR: expected $THR but it is missing/empty — see the GISTIC log above." >&2
  exit 1
fi
