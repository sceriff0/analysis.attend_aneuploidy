#!/usr/bin/env bash
# =============================================================================
# fetch_tcga_ucec.sh — download the TCGA UCEC clinical table for report 12
#
# Cohort B (PRIMARY tumours) for the ATTEND-vs-TCGA MMRd aneuploidy contrast. Pulls the
# cBioPortal PanCancer-Atlas study `ucec_tcga_pan_can_atlas_2018` and extracts ONLY the
# two clinical files the loader needs (sample + patient), which carry the per-sample
# ANEUPLOIDY_SCORE (Taylor et al. 2018), SUBTYPE (UCEC_MSI/_POLE/_CN_HIGH/_CN_LOW),
# MSI_SCORE_MANTIS, FRACTION_GENOME_ALTERED, TMB, and OS survival. load_tcga_ucec()
# reads them; the full tarball (~200 MB) is not kept.
#
# For TP53 status in the reproduction check, also extract data_mutations.txt (large) —
# uncomment the line below.
# =============================================================================
set -euo pipefail

DEST="${DEST:-data/tcga}"
STUDY="ucec_tcga_pan_can_atlas_2018"
URL="https://cbioportal-datahub.s3.amazonaws.com/${STUDY}.tar.gz"

mkdir -p "$DEST"
echo "Downloading $STUDY clinical files from cBioPortal datahub..."
# stream the tarball and extract only the clinical files (and optionally mutations)
curl -fsSL "$URL" | tar -xz -C "$DEST" --strip-components=1 \
  "${STUDY}/data_clinical_sample.txt" \
  "${STUDY}/data_clinical_patient.txt"
  # "${STUDY}/data_mutations.txt"   # <- uncomment for TP53 status (large)

echo "Wrote:"
ls -1 "$DEST"/data_clinical_*.txt
echo ""
echo "Now in R:  source('code/load_wes_results.R'); load_tcga_ucec()  # -> harmonised TCGA table"
echo "Then knit analysis/12_mmrd_aneuploidy_crosscohort.Rmd"

# --- Reproducible R alternative (no shell / different mirror) -----------------
# If you prefer Bioconductor:
#   BiocManager::install("cBioPortalData")
#   mae <- cBioPortalData::cBioDataPack("ucec_tcga_pan_can_atlas_2018")
#   cd  <- as.data.frame(MultiAssayExperiment::colData(mae))
#   readr::write_tsv(cd, "data/tcga/data_clinical_sample.txt")   # colData already merged
# (the loader's #-comment skip is a no-op on this flat file — it still reads.)
