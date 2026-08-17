# =============================================================================
# load_wes_results.R  —  WES-DERIVED GENOMIC LOADERS
#
# Sourced by report 01, 31, 33, 34, 35, 40 and 41. All tables are keyed in BARCODE space and
# mapped to pid via the barcode crosswalk. The barcode-suffix stripping below
# (`-1TAD104`, `_tumor_only`) normalises the sequencing IDs back to the bare
# barcode that the gianlu crosswalk (TUMOR_BARCODE) uses.
# Raw inputs live in data/ (synced from HPC: /hpcnfs/scratch/P_DIMA_ATTEND).
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(fs)
  library(tidyverse)
  library(data.table)
})

# aneuploidy: data/aneu_scores.csv — returned as-is (must carry `sequencing_id`).
load_aneuploidy_scores <- function() {
  read.csv(here("data", "aneu_scores.csv")) |> as_tibble()
}

# HRD: data/hrd_scores.csv — strip the sequencing suffix so `Sample` is the bare barcode.
load_hrd_scores <- function() {
  read.csv(here("data", "hrd_scores.csv")) |>
    as_tibble() |>
    mutate(Sample = str_remove(Sample, "-1TAD104_tumor_only|_tumor_only"))
}

# MSI: data/msi_scores.csv — one JSON-ish line per sample; pull barcode + % unstable.
load_msi_scores <- function() {
  tibble(line = read_lines(here("data", "msi_scores.csv"))) |>
    filter(str_detect(line, "_tumor_only")) |>
    extract(
      line,
      into   = c("ID", "MSI_SCORE"),
      regex  = '^(.+?)_+tumor_only.*?PercentageUnstableSites"?\\s*:\\s*"?([0-9.]+)',
      remove = TRUE
    ) |>
    mutate(ID        = str_remove(ID, "-1TAD104$"),
           MSI_SCORE = as.numeric(MSI_SCORE)) |>
    select(ID, MSI_SCORE)
}

# TMB: data/tmb_scores.csv — one line per sample; pull barcode + TMB (mut/Mb).
load_tmb_scores <- function() {
  tibble(line = read_lines(here("data", "tmb_scores.csv"))) |>
    filter(str_detect(line, "_tumor_only")) |>
    extract(
      line,
      into   = c("ID", "TMB_SCORE"),
      regex  = '^(.+?)_+tumor_only.*?TMB,([0-9.]+)',
      remove = TRUE
    ) |>
    mutate(ID        = str_remove(ID, "-1TAD104$"),
           TMB_SCORE = as.numeric(TMB_SCORE)) |>
    select(ID, TMB_SCORE)
}

# --- MAF -> per-sample gene status ------------------------------------------
# Each .maf (one per sample, ~399 annotation columns) is reduced to one row:
# ID + <GENE>_status, where a gene is "ABNORMAL" if any of its variant rows is
# called pathogenic. Consumed by pathogenic_by_patient() in attend_classes.R
# (never mean-collapsed). Add genes via `genes` — columns follow <GENE>_status.
#
# Pathogenicity is read from the AUTHORITATIVE ClinVar significance column only
# (MAF_SIG_COLS). This is deliberate: scanning all columns over-calls, because
# unrelated fields contain the substring (e.g. `Pathway`, `clinvar: Clinvar`) and
# `ClinVar_VCF_CLNSIGCONF` holds "Conflicting_classifications_of_pathogenicity",
# which is NOT a pathogenic call. The (?![a-z]) lookahead excludes that "...icity".
#
# To change the definition (these are SOMATIC MAFs, so ClinVar alone can under-call
# novel truncating variants): add "am_class" (AlphaMissense) to MAF_SIG_COLS, or
# switch to Variant_Classification-based calling against attend_pathogenic_classes.
genes          <- c("TP53")
MAF_SIG_COLS   <- c("ClinVar_VCF_CLNSIG")
MAF_PATHOGENIC <- regex("(likely_)?pathogenic(?![a-z])", ignore_case = TRUE)

gene_status <- function(maf, gene, sig_cols = MAF_SIG_COLS, pattern = MAF_PATHOGENIC) {
  rows <- maf |> filter(Hugo_Symbol == gene)
  if (nrow(rows) == 0) return("NORMAL")
  cols <- intersect(sig_cols, names(rows))
  if (length(cols) == 0) return(NA_character_)              # annotation column absent
  hit <- rows |> filter(if_any(all_of(cols), ~ str_detect(as.character(.x), pattern)))
  if (nrow(hit) > 0) "ABNORMAL" else "NORMAL"
}

process_maf <- function(maf_path) {
  id  <- path_ext_remove(path_file(maf_path))
  maf <- fread(maf_path, skip = "Hugo_Symbol", sep = "\t", header = TRUE,
               quote = "", fill = Inf, na.strings = c(".", "", "NA")) |>
    as_tibble()

  status <- set_names(map_chr(genes, ~ gene_status(maf, .x)), paste0(genes, "_status"))
  tibble(ID = str_remove(id, "-1TAD104"), !!!status)
}

# Combine every per-sample MAF into the wide table (ID + <GENE>_status columns).
# On the cluster the .maf files live in `variant_annotation`; sync them into
# data/variant_annotations/ (or repoint maf_folder) before running.
load_maf_data <- function() {
  maf_folder <- here("data", "variant_annotations")
  # All files in the folder are treated as MAFs (original behaviour). If the dir
  # mixes in non-MAF files, add a glob (e.g. glob = "*.maf") to dir_ls().
  dir_ls(maf_folder) |>
    map(process_maf) |>
    list_rbind()
}

# --- Long per-variant MAF (for maftools / oncoplots) ------------------------
# Reads every annotated MAF and returns the standard columns maftools needs, with
# Tumor_Sample_Barcode set from the file name (stripped) so it maps to the barcode
# crosswalk. One row per variant across all samples. Reads all 247 files — slow.
maf_standard_cols <- c(
  "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
  "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type"
)

read_one_maf_long <- function(maf_path) {
  id <- str_remove(path_ext_remove(path_file(maf_path)), "-1TAD104")
  fread(maf_path, skip = "Hugo_Symbol", sep = "\t", header = TRUE, quote = "",
        fill = Inf, na.strings = c(".", "", "NA")) |>
    as_tibble() |>
    select(any_of(maf_standard_cols)) |>
    mutate(Tumor_Sample_Barcode = id)
}

load_maf_long <- function() {
  here("data", "variant_annotations") |>
    dir_ls() |>
    map(read_one_maf_long) |>
    list_rbind()
}

# --- Long MAF for ANCESTRY-AWARE TMB (report 01) ----------------------------
# Keeps the columns tmb_variant_table() needs: consequence, protein change (for
# definitional POLE hotspots), and gnomAD exome v4.1 per-sub-population + grpmax
# allele frequencies. Tumor_Sample_Barcode is set from the file name (stripped)
# so it maps via the barcode crosswalk. Returns empty tibble if the folder is
# absent — knit-safe.
maf_tmb_cols <- c(
  "Hugo_Symbol", "Variant_Classification",
  "gnomAD_exome_AF_afr", "gnomAD_exome_AF_amr", "gnomAD_exome_AF_asj",
  "gnomAD_exome_AF_eas", "gnomAD_exome_AF_fin", "gnomAD_exome_AF_nfe",
  "gnomAD_exome_AF_sas", "gnomAD_exome_AF_mid", "gnomAD_exome_AF_grpmax"
)

load_maf_tmb <- function() {
  folder <- here("data", "variant_annotations")
  if (!dir_exists(folder)) {
    message("load_maf_tmb(): no folder ", folder, " — returning empty."); return(tibble())
  }
  dir_ls(folder) |>
    map(function(p) {
      id <- str_remove(path_ext_remove(path_file(p)), "-1TAD104")
      fread(p, skip = "Hugo_Symbol", sep = "\t", header = TRUE, quote = "",
            fill = Inf, na.strings = c(".", "", "NA")) |>
        as_tibble() |>
        select(any_of(maf_tmb_cols)) |>
        mutate(across(starts_with("gnomAD_exome_AF"), ~ suppressWarnings(as.numeric(.x))),
               Tumor_Sample_Barcode = id)
    }) |>
    list_rbind()
}

# Derive a ClassifyCNV sample id from its file path — the SAME transform load_cnv_data()
# applies inside its map() body, factored out so the two can never drift. Exists because
# `cc$id_strip` (e.g. "-1TAD104|_tumor_only") strips MULTIPLE sequencing-id suffixes down
# to one bare barcode BY DESIGN (the same sample can be resequenced under a different
# suffix) — so `length(files)` overcounts whenever that happens, and `dir_ls(recurse=TRUE)`
# overcounts further if the same sample's file also appears under two subfolders. Distinct
# ids, not file count, is what "how many samples were profiled" means.
.classifycnv_file_id <- function(p, cc) {
  id <- str_remove(path_file(p), "\\.cnv\\.annotated\\.tsv$")
  str_remove(id, cc$id_strip)
}

# Count of samples PROFILED among a set of ClassifyCNV file paths — i.e. distinct ids
# derived by .classifycnv_file_id(), not length(files) (see that function's comment for
# why the two differ: id_strip collapses a resequenced sample's two files onto one id).
# Factored out of load_cnv_data() so it is testable on its own: this is pure string
# manipulation plus n_distinct(), no I/O, so it runs without `fs`/`data.table`/`here` —
# unlike load_cnv_data() itself, which this repro's test suite cannot execute locally
# (see test_classifycnv_file_id.R). load_cnv_data() calls this for attr(out, "n_profiled").
.classifycnv_n_profiled <- function(files, cc) {
  dplyr::n_distinct(vapply(files, .classifycnv_file_id, character(1), cc = cc))
}

# --- CNV: ClassifyCNV-annotated segment calls (report 06) -------------------
# variantalker output: one file per sample, one row per ALTERED segment (DUP/DEL)
# with `logRatio` (log2 ratio) + `CNF` (= 2^logRatio). The ACMG `1A..5H` scoring
# columns are clinical-pathogenicity criteria, irrelevant to broad-SCNA clustering,
# so we keep only location + Type + logRatio/CNF. ID comes from the file name
# (sequencing suffix stripped) so it maps to the barcode crosswalk.
# Sync the HPC tree into data/cnv_annotations/<s>/<s>.cnv.annotated.tsv first.
# Returns an EMPTY tibble (not an error) when the folder is absent — knit-safe.
load_cnv_data <- function(cnv = attend_cnv) {
  cc     <- cnv$classifycnv
  folder <- here("data", cc$dir)
  if (!dir_exists(folder)) {
    message("load_cnv_data(): no folder ", folder, " — returning empty."); return(tibble())
  }
  files <- dir_ls(folder, recurse = TRUE, type = "file", glob = "*.cnv.annotated.tsv")
  if (length(files) == 0) {
    message("load_cnv_data(): no *.cnv.annotated.tsv under ", folder)
    out <- tibble()
    attr(out, "n_profiled") <- 0L
    return(out)
  }
  out <- files |>
    map(function(p) {
      id <- .classifycnv_file_id(p, cc)
      dt <- fread(p, sep = "\t", header = TRUE, quote = "",
                  na.strings = c(".", "", "NA")) |> as_tibble()
      keep <- intersect(c(cc$chrom_col, cc$start_col, cc$end_col,
                          cc$type_col, cc$value_col, cc$cnf_col), names(dt))
      dt |>
        select(all_of(keep)) |>
        rename(any_of(c(Chromosome = cc$chrom_col, Start = cc$start_col,
                        End = cc$end_col, Type = cc$type_col,
                        logRatio = cc$value_col, CNF = cc$cnf_col))) |>
        mutate(ID = id)
    }) |>
    list_rbind() |>
    relocate(ID)
  # Attach AFTER the final dplyr verb — relocate() (like every dplyr verb) drops
  # non-standard attributes, so setting this any earlier would lose it. .classifycnv_n_profiled()
  # counts DISTINCT derived ids (via the same .classifycnv_file_id() the map() body used
  # above), not length(files) — two files can legitimately derive the same id (see that
  # helper's comment), and length(files) would then overcount samples PROFILED, inflating
  # the segment_pileup() denominator and deflating every reported frequency.
  attr(out, "n_profiled") <- .classifycnv_n_profiled(files, cc)
  out
}

# --- CNV: DRAGEN genome-wide segmentation (.seg) for the Fig-1a heatmap ------
# One .seg per sample under data/seg/ (subfolders OK) — DRAGEN's genome-wide segmented
# copy number (the same file that feeds ASCETS), the input to the TCGA Fig-1a-style
# GENOME-BINNED clustering (segments_to_bin_matrix() in attend_classes.R). Returns long
# segments: ID + Chromosome + Start + End + Segment_Mean, with ID taken from the file
# name (suffix stripped) so it maps via the barcode crosswalk. Column names are
# AUTO-DETECTED across DRAGEN variants (attend_cnv$seg$*_col); a file whose columns
# can't be resolved is skipped with a message. Empty tibble when data/seg/ is absent
# (knit-safe). NOTE: non-ATTEND-barcode files (controls, other cohorts) load fine but
# are dropped later by the barcode crosswalk since they have no pid.
load_seg_data <- function(cnv = attend_cnv) {
  sc     <- cnv$seg
  folder <- here("data", sc$dir)
  if (!dir_exists(folder)) {
    message("load_seg_data(): no folder ", folder, " — returning empty."); return(tibble())
  }
  files <- dir_ls(folder, recurse = TRUE, type = "file", glob = "*.seg")
  if (length(files) == 0) {
    message("load_seg_data(): no *.seg under ", folder); return(tibble())
  }
  pick <- function(nms, cands) { hit <- intersect(cands, nms); if (length(hit)) hit[[1]] else NA_character_ }
  files |>
    map(function(p) {
      id <- str_remove(path_ext_remove(path_file(p)), sc$id_strip)
      dt <- tryCatch(
        fread(p, header = TRUE, na.strings = c(".", "", "NA"), fill = TRUE) |> as_tibble(),
        error = function(e) { message("load_seg_data(): unreadable ", path_file(p), " — ", conditionMessage(e)); NULL })
      if (is.null(dt) || nrow(dt) == 0) return(NULL)
      cc <- pick(names(dt), sc$chrom_col); ss <- pick(names(dt), sc$start_col)
      ee <- pick(names(dt), sc$end_col);   vv <- pick(names(dt), sc$value_col)
      if (any(is.na(c(cc, ss, ee, vv)))) {
        message("load_seg_data(): could not resolve columns in ", path_file(p),
                " (have: ", paste(names(dt), collapse = ", "), ") — skipping."); return(NULL)
      }
      tibble(ID           = id,
             Chromosome   = as.character(dt[[cc]]),
             Start        = suppressWarnings(as.numeric(dt[[ss]])),
             End          = suppressWarnings(as.numeric(dt[[ee]])),
             Segment_Mean = suppressWarnings(as.numeric(dt[[vv]])))
    }) |>
    purrr::compact() |>
    list_rbind()
}

# --- TCGA UCEC reference cohort (legacy PanCanAtlas export; fallback) ---------
# Reads the cBioPortal UCEC PanCancer-Atlas clinical export from data/tcga/ and returns a
# harmonised per-sample tibble for the ATTEND(metastatic)-vs-TCGA(primary) contrast:
#   pid, aneuploidy, subtype, mmrd (logical), fga, tmb, os_months, os_event (1/0)
# The cBioPortal data_clinical_{sample,patient}.txt files have 4 leading #-comment lines
# before the real header; we skip them by finding the first non-# line. Returns an empty
# tibble when data/tcga/ is absent — knit-safe. Get the files via code/fetch_tcga_ucec.R.
read_cbioportal_clinical <- function(path) {
  ln <- readr::read_lines(path)
  hdr <- which(!startsWith(ln, "#"))[1]          # first non-comment line = real header
  readr::read_tsv(I(paste(ln[hdr:length(ln)], collapse = "\n")), show_col_types = FALSE)
}

load_tcga_ucec <- function(ref = attend_tcga_ref) {
  folder <- here("data", ref$dir)
  if (!dir_exists(folder)) { message("load_tcga_ucec(): no folder ", folder, " — returning empty."); return(tibble()) }
  fs_ <- dir_ls(folder, recurse = TRUE, type = "file", glob = ref$sample_glob)
  fp_ <- dir_ls(folder, recurse = TRUE, type = "file", glob = ref$patient_glob)
  if (length(fs_) == 0) { message("load_tcga_ucec(): no ", ref$sample_glob, " under ", folder); return(tibble()) }

  samp <- tryCatch(read_cbioportal_clinical(fs_[[1]]), error = function(e) NULL)
  if (is.null(samp)) { message("load_tcga_ucec(): could not read ", fs_[[1]]); return(tibble()) }
  if (length(fp_) > 0) {
    pat <- tryCatch(read_cbioportal_clinical(fp_[[1]]), error = function(e) NULL)
    if (!is.null(pat) && ref$id_col %in% names(samp) && ref$id_col %in% names(pat))
      samp <- dplyr::left_join(samp, pat, by = ref$id_col, suffix = c("", ".pat"))
  }

  g <- function(col) if (col %in% names(samp)) samp[[col]] else rep(NA, nrow(samp))
  subtype <- as.character(g(ref$subtype_col))
  mantis  <- suppressWarnings(as.numeric(g(ref$msi_score_col)))
  mmrd    <- (!is.na(subtype) & str_detect(subtype, ref$msi_subtype_pattern)) |
             (!is.na(mantis)  & mantis > ref$msi_hi_mantis)
  osraw   <- as.character(g(ref$os_event_col))
  tibble(
    pid        = as.character(g(ref$id_col)),
    aneuploidy = suppressWarnings(as.numeric(g(ref$aneuploidy_col))),
    subtype    = subtype,
    mmrd       = mmrd,
    fga        = suppressWarnings(as.numeric(g(ref$fga_col))),
    tmb        = suppressWarnings(as.numeric(g(ref$tmb_col))),
    os_months  = suppressWarnings(as.numeric(g(ref$os_time_col))),
    os_event   = dplyr::case_when(str_detect(osraw, "^1|DECEASED") ~ 1L,
                                  str_detect(osraw, "^0|LIVING")   ~ 0L,
                                  TRUE ~ NA_integer_)
  ) |> dplyr::filter(!is.na(pid))
}

# --- TCGA UCEC 2013 publication cohort (report 07, TMB reference) -------------
# Reads data/tcga_2013/ (populated by code/fetch_tcga_ucec_2013.R) and returns ONE ROW PER
# SAMPLE — because that is the grain the Fig-1a heatmap needs: its columns are tumours
# (sample barcodes, the key of the .seg file), while the paper's SUBTYPE call is recorded
# per PATIENT. This loader resolves that mismatch in one place by broadcasting the
# patient-level columns (SUBTYPE, survival) down onto each of that patient's samples, so
# every downstream join is a plain join on `ID`.
#
# Columns: ID (sample barcode) / pid (patient) / subtype (ordered factor, cascade order) /
#   subtype_short / cn_cluster_k4 (PUBLISHED Fig-1a cluster, factor 1-4) / mmrd (logical) /
#   msi_call / mlh1_silencing / mutrate_cluster / fga / tmb / aneuploidy (retrofitted from
#   PanCancer-Atlas, NA where absent) / os_months / os_event / pfs_months / pfs_event.
#
# `mmrd` follows the same two-source convention as load_tcga_ucec(): SUBTYPE contains "MSI"
# OR the study's own 7-marker call is MSI-H. Empty tibble when the folder is absent — knit-safe.
load_tcga_ucec_2013 <- function(ref = attend_tcga_ref_2013) {
  folder <- here("data", ref$dir)
  if (!dir_exists(folder)) {
    message("load_tcga_ucec_2013(): no folder ", folder,
            " — run `Rscript code/fetch_tcga_ucec_2013.R`. Returning empty."); return(tibble())
  }
  fs_ <- dir_ls(folder, recurse = TRUE, type = "file", glob = ref$sample_glob)
  fp_ <- dir_ls(folder, recurse = TRUE, type = "file", glob = ref$patient_glob)
  if (length(fs_) == 0) { message("load_tcga_ucec_2013(): no ", ref$sample_glob, " under ", folder); return(tibble()) }

  samp <- tryCatch(read_cbioportal_clinical(fs_[[1]]), error = function(e) NULL)
  if (is.null(samp)) { message("load_tcga_ucec_2013(): could not read ", fs_[[1]]); return(tibble()) }
  if (length(fp_) > 0) {
    pat <- tryCatch(read_cbioportal_clinical(fp_[[1]]), error = function(e) NULL)
    if (!is.null(pat) && ref$id_col %in% names(samp) && ref$id_col %in% names(pat))
      samp <- dplyr::left_join(samp, pat, by = ref$id_col, suffix = c("", ".pat"))
  }
  # Aneuploidy is retrofitted from the PanCancer-Atlas study, joined on PATIENT_ID.
  fa <- file.path(folder, ref$aneuploidy_file)
  if (file.exists(fa)) {
    aneu <- tryCatch(readr::read_tsv(fa, show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(aneu) && all(c(ref$id_col, ref$aneuploidy_col) %in% names(aneu)))
      samp <- dplyr::left_join(samp, aneu[, c(ref$id_col, ref$aneuploidy_col)], by = ref$id_col)
  } else message("load_tcga_ucec_2013(): no ", ref$aneuploidy_file,
                 " — aneuploidy will be all-NA (re-run the fetcher).")

  g   <- function(col) if (col %in% names(samp)) samp[[col]] else rep(NA, nrow(samp))
  ev  <- function(col, pos) {                      # "1:DECEASED" / "0:LIVING" -> 1/0
    raw <- as.character(g(col))
    dplyr::case_when(str_detect(raw, paste0("^1|", pos)) ~ 1L,
                     str_detect(raw, "^0")               ~ 0L,
                     TRUE ~ NA_integer_)
  }
  subtype <- as.character(g(ref$subtype_col))
  msi     <- as.character(g(ref$msi_call_col))
  # --- MMRd -----------------------------------------------------------------
  # Two rules, selected by ref$mmrd_rule. The DEFAULT is "tcga_exact".
  #
  # "tcga_exact" — MMRd IS the paper's MSI (Hyper-mutated) integrative cluster and nothing
  #   else (n = 65). Every other CLASSIFIED tumour (POLE / CN-high / CN-low) is MMRp, and the
  #   141 tumours TCGA never classified are NA — unclassified, not proficient. This is the
  #   only rule whose output can be quoted as "TCGA's MSI group".
  #
  # "broad" — additionally accepts a 7-marker MSI-H call in tumours with no integrative
  #   cluster (n = 124, i.e. +59). Broader coverage, but NOT TCGA's classification.
  #
  # ⚠️ POLE TAKES PRECEDENCE under either rule. The cascade is POLE -> MSI -> CN-low/CN-high,
  # so a POLE-ultramutated tumour is never MMRd even though 3 of them are MSI-H by the marker
  # assay (TCGA-AP-A051, TCGA-AP-A059, TCGA-D1-A103) — ultramutation destabilises
  # microsatellites as a consequence of mutation burden, not through MMR loss. Under
  # "tcga_exact" the SUBTYPE column already encodes this; the guard matters for "broad".
  pole_val <- !is.na(subtype) & str_detect(subtype, ref$pole_subtype_pattern)
  sub_msi  <- !is.na(subtype) & str_detect(subtype, ref$msi_subtype_pattern)
  if (identical(ref$mmrd_rule, "broad")) {
    mmrd_pos <- !pole_val & (sub_msi | (!is.na(msi) & msi == "MSI-H"))
    mmrd_inf <- !is.na(subtype) | (!is.na(msi) & msi %in% c("MSS", "MSI-L"))
    mmrd_val <- ifelse(mmrd_pos, TRUE, ifelse(mmrd_inf, FALSE, NA))
  } else {
    # TCGA-exact: the integrative cluster is the sole authority. Classified -> TRUE/FALSE;
    # unclassified -> NA. The marker call is deliberately NOT consulted.
    mmrd_val <- ifelse(is.na(subtype), NA, sub_msi)
  }

  tibble(
    ID              = as.character(g(ref$sample_id_col)),
    pid             = as.character(g(ref$id_col)),
    subtype         = factor(subtype, levels = ref$subtype_levels),
    subtype_short   = factor(unname(ref$subtype_short[subtype]),
                             levels = unname(ref$subtype_short)),
    cn_cluster_k4   = factor(as.character(g(ref$cluster_col)), levels = as.character(1:4)),
    # TRUE for TCGA's 232-sample multiplatform core set. The integrated SUBTYPE exists for
    # exactly these samples and for no others, so this is the honest way to restrict to the
    # classified cohort — see the comment on core_col in attend_classes.R.
    core_sample     = toupper(trimws(as.character(g(ref$core_col)))) %in% c("Y", "YES", "TRUE", "1"),
    mmrd            = mmrd_val,
    pole            = pole_val,          # POLE-ultramutated; excluded from mmrd by the cascade
    msi_call        = msi,
    # Grade and stage are complete in this cohort (373/373) and are the two clinical
    # covariates the aneuploidy panels colour by. Ordered factors so the legend and any
    # model term run low->high rather than alphabetically. "Unknown" stage (n = 3) becomes
    # NA rather than a 5th level — it is missingness, not a stage between IV and I.
    grade           = factor(as.character(g(ref$grade_col)),
                             levels = ref$grade_levels, ordered = TRUE),
    stage           = factor(ifelse(as.character(g(ref$stage_col)) == ref$stage_unknown,
                                    NA_character_, as.character(g(ref$stage_col))),
                             levels = ref$stage_levels, ordered = TRUE),
    mlh1_silencing  = suppressWarnings(as.integer(g(ref$mlh1_col))),
    mutrate_cluster = as.character(g(ref$mutrate_col)),
    fga             = suppressWarnings(as.numeric(g(ref$fga_col))),
    tmb             = suppressWarnings(as.numeric(g(ref$tmb_col))),
    aneuploidy      = suppressWarnings(as.numeric(g(ref$aneuploidy_col))),
    os_months       = suppressWarnings(as.numeric(g(ref$os_time_col))),
    os_event        = ev(ref$os_event_col,  "DECEASED"),
    pfs_months      = suppressWarnings(as.numeric(g(ref$dfs_time_col))),
    pfs_event       = ev(ref$dfs_event_col, "Recurred|Progressed")
  ) |> dplyr::filter(!is.na(ID))
}

# Long segment table for the 2013 cohort — the SAME shape load_seg_data() returns for
# ATTEND (ID / Chromosome / Start / End / Segment_Mean), so segments_to_bin_matrix() and
# segments_to_arm_matrix() consume it unchanged. ID is the SAMPLE barcode, matching
# load_tcga_ucec_2013()$ID. Coordinates are hg19 — pass load_arm_boundaries(build = "hg19").
# Empty tibble when the folder/file is absent — knit-safe.
load_tcga_2013_seg <- function(ref = attend_tcga_ref_2013) {
  folder <- here("data", ref$dir)
  if (!dir_exists(folder)) { message("load_tcga_2013_seg(): no folder ", folder); return(tibble()) }
  files <- dir_ls(folder, recurse = TRUE, type = "file", glob = ref$seg_glob)
  if (length(files) == 0) {
    message("load_tcga_2013_seg(): no ", ref$seg_glob, " under ", folder,
            " — run `Rscript code/fetch_tcga_ucec_2013.R`."); return(tibble())
  }
  dt <- tryCatch(fread(files[[1]], header = TRUE, na.strings = c(".", "", "NA")) |> as_tibble(),
                 error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) { message("load_tcga_2013_seg(): unreadable ", files[[1]]); return(tibble()) }
  need <- c("Sample", "Chromosome", "Start", "End", "Segment_Mean")
  if (!all(need %in% names(dt))) {
    message("load_tcga_2013_seg(): expected ", paste(need, collapse = "/"),
            " but have ", paste(names(dt), collapse = ", ")); return(tibble())
  }
  dt |> transmute(ID           = as.character(Sample),
                  # defensive: the fetcher already maps 23/24 -> X/Y, but a hand-placed
                  # .seg from another source may not have been normalised.
                  Chromosome   = dplyr::recode(str_remove(as.character(Chromosome), "^chr"),
                                               "23" = "X", "24" = "Y"),
                  Start        = suppressWarnings(as.numeric(Start)),
                  End          = suppressWarnings(as.numeric(End)),
                  Segment_Mean = suppressWarnings(as.numeric(Segment_Mean)))
}

# ASCETS run on the TCGA 2013 segments — the aneuploidy score that is actually COMMENSURATE
# with ATTEND's. Produced by code/run_ascets_tcga.sh; reads data/tcga_2013/ascets/.
#
# WHY THIS AND NOT THE cBIOPORTAL COLUMN. cBioPortal's ANEUPLOIDY_SCORE is Taylor et al.
# (Cancer Cell 2018): a COUNT of altered arms out of 39 autosomal arms (integer, 0-39).
# ATTEND's is ASCETS: the FRACTION of *evaluable* arms altered,
#   #{call in AMP,DEL} / #{call != LOWCOV}   (0-1)
# — a different denominator, and a per-sample one, so the two cannot be placed on a common
# axis by rescaling. Running ASCETS on TCGA's own segments produces the same statistic for
# both cohorts. Report 16 prefers this whenever it exists and falls back to the cBioPortal
# column otherwise.
#
# Returns ID (sample barcode) + ascets_aneuploidy, plus — when the arm-level calls file is
# present — the diagnostic columns that decide whether the comparison is trustworthy:
#   n_arms_evaluable  arms with a call other than LOWCOV (the score's DENOMINATOR)
#   n_arms_altered    arms called AMP or DEL             (the score's NUMERATOR)
#   n_arms_lowcov     arms dropped for insufficient coverage
# A cohort whose n_arms_evaluable varies widely has a moving denominator, which inflates
# the score for poorly covered samples — surface it rather than let it pass silently.
# Empty tibble when the folder is absent — knit-safe.
load_tcga_2013_ascets <- function(ref = attend_tcga_ref_2013) {
  ac     <- ref$ascets
  folder <- here("data", ac$dir)
  if (!dir_exists(folder)) {
    message("load_tcga_2013_ascets(): no folder ", folder,
            " — run `bash code/run_ascets_tcga.sh`. Returning empty."); return(tibble())
  }
  fs_ <- dir_ls(folder, recurse = TRUE, type = "file", glob = ac$aneuploidy_glob)
  if (length(fs_) == 0) {
    message("load_tcga_2013_ascets(): no ", ac$aneuploidy_glob, " under ", folder); return(tibble())
  }
  sc <- tryCatch(fread(fs_[[1]], sep = "\t", header = TRUE) |> as_tibble(), error = function(e) NULL)
  if (is.null(sc) || nrow(sc) == 0) { message("load_tcga_2013_ascets(): unreadable ", fs_[[1]]); return(tibble()) }
  names(sc)[1] <- "ID"
  num <- names(sc)[map_lgl(sc, is.numeric)]
  if (length(num) == 0) { message("load_tcga_2013_ascets(): no numeric score column."); return(tibble()) }
  out <- sc |> transmute(ID = as.character(ID),
                         ascets_aneuploidy = suppressWarnings(as.numeric(.data[[num[1]]])))

  # arm-level calls -> the numerator/denominator behind each score
  fc_ <- dir_ls(folder, recurse = TRUE, type = "file", glob = ac$calls_glob)
  if (length(fc_) > 0) {
    cl <- tryCatch(fread(fc_[[1]], sep = "\t", header = TRUE, na.strings = c("", "NA")) |> as_tibble(),
                   error = function(e) NULL)
    if (!is.null(cl) && nrow(cl) > 0) {
      names(cl)[1] <- "ID"
      diag <- cl |>
        pivot_longer(-ID, names_to = "arm", values_to = "call") |>
        mutate(call = as.character(call)) |>
        group_by(ID = as.character(ID)) |>
        summarise(n_arms_evaluable = sum(!is.na(call) & call != "LOWCOV"),
                  n_arms_altered   = sum(!is.na(call) & call %in% c("AMP", "DEL")),
                  n_arms_lowcov    = sum(!is.na(call) & call == "LOWCOV"),
                  .groups = "drop")
      out <- dplyr::left_join(out, diag, by = "ID")
    }
  }
  out
}

# --- TCGA's PUBLISHED GISTIC peaks -------------------------------------------
# NOTE: currently UNUSED by any report — its caller was removed with the
# published-peak clustering report. Kept: the peak file is bundled and the loader is
# tested, so re-wiring it is a one-line change.
# Reads data/tcga_gistic_peaks_2013.tsv (Suppl. Data File S2.1 of nature12113, "all tumor"
# sheets, converted by hand once). Long format: one row per peak x gene.
#   peak_id / direction (amp|del) / cytoband / q_value / residual_q / chrom / start / end / gene
# Coordinates are hg19 — the build of the TCGA segments, NOT of ATTEND's.
# Returns an empty tibble when the file is absent — knit-safe.
load_tcga_published_peaks <- function(ref = attend_tcga_ref_2013) {
  p <- here("data", ref$published_peaks)
  if (!file.exists(p)) {
    message("load_tcga_published_peaks(): missing ", p); return(tibble())
  }
  readr::read_tsv(p, show_col_types = FALSE, col_types = readr::cols(.default = "c")) |>
    mutate(start = suppressWarnings(as.numeric(start)),
           end   = suppressWarnings(as.numeric(end)),
           q_value = suppressWarnings(as.numeric(q_value)))
}

# Sample x PEAK matrix for a cohort whose segments are on the SAME BUILD as the peaks
# (hg19 -> the TCGA 2013 .seg). For each sample and each published peak, the length-weighted
# mean log2 seg-mean across the peak's wide boundaries; peaks with no overlapping segment are
# 0 (copy-neutral). This is the COORDINATE route: it needs no GISTIC run at all, because the
# peaks are already given. `seg` is the long table from load_tcga_2013_seg().
#
# `threshold`: NULL keeps the continuous length-weighted log2. A numeric c(low, high)
# DISCRETISES to {-2,-1,0,+1,+2} — |x| < low -> 0, low <= |x| < high -> +-1, |x| >= high -> +-2.
# The paper's Suppl. Methods say clustering was run on "THRESHOLDED relative copy number data
# in significantly reoccurring ... regions", so the discrete form is the faithful one; it also
# matches ATTEND's route, whose GISTIC gene calls are already thresholded. The paper does not
# state its cut points, hence the parameter.
published_peak_matrix_from_seg <- function(seg, peaks, threshold = NULL) {
  if (is.null(seg) || nrow(seg) == 0 || is.null(peaks) || nrow(peaks) == 0) return(NULL)
  pk <- peaks |>
    distinct(peak_id, chrom, start, end) |>
    filter(!is.na(start), !is.na(end), end > start)
  if (!nrow(pk)) return(NULL)
  sg <- seg |>
    transmute(ID = as.character(ID), chrom = str_remove(as.character(Chromosome), "^chr"),
              s = as.numeric(Start), e = as.numeric(End), v = as.numeric(Segment_Mean)) |>
    filter(!is.na(ID), !is.na(s), !is.na(e), !is.na(v), e > s)
  pk |>
    inner_join(sg, by = "chrom", relationship = "many-to-many") |>
    mutate(ov = pmax(0, pmin(e, end) - pmax(s, start))) |>
    filter(ov > 0) |>
    group_by(ID, peak_id, width = end - start) |>
    summarise(v = sum(ov * v) / first(width), .groups = "drop") |>
    select(ID, peak_id, v) |>
    tidyr::complete(ID, peak_id = pk$peak_id, fill = list(v = 0)) |>
    mutate(v = if (is.null(threshold)) v else {
      lo <- min(threshold); hi <- max(threshold)
      sign(v) * (as.integer(abs(v) >= lo) + as.integer(abs(v) >= hi))
    }) |>
    pivot_wider(names_from = peak_id, values_from = v, values_fill = 0)
}

# Sample x GENE matrix restricted to the published peak genes — the GENE route, for a cohort
# on a DIFFERENT build (ATTEND is hg38). Gene symbols are build-independent, so no liftover
# is needed: we simply subset that cohort's own GISTIC all_thresholded.by_genes.txt to the
# genes TCGA reported inside its peaks. `thr` is the tibble from load_gistic_thresholded()
# (ID + one column per gene, values in {-2..+2}); it may already be peak-restricted by that
# cohort's OWN peaks — this replaces that restriction with TCGA's.
# `collapse = TRUE` (default) then averages the genes WITHIN each peak, returning a sample x
# PEAK matrix on the same 79 columns as published_peak_matrix_from_seg() — so the two cohorts
# occupy one shared feature space and their clusterings are directly comparable. Set FALSE to
# keep gene-level columns. NOTE the values still differ in KIND between the routes: thresholded
# GISTIC calls averaged per peak here, length-weighted mean log2 there. Comparable in structure
# and sign, not in units — z-score within cohort before comparing feature values across them.
published_peak_matrix_from_genes <- function(thr, peaks, id_col = "ID", collapse = TRUE) {
  if (is.null(thr) || nrow(thr) == 0 || is.null(peaks) || nrow(peaks) == 0) return(NULL)
  genes <- intersect(setdiff(names(thr), id_col), unique(peaks$gene))
  if (!length(genes)) {
    message("published_peak_matrix_from_genes(): no overlap between the cohort's GISTIC ",
            "genes and the ", dplyr::n_distinct(peaks$gene), " published peak genes — ",
            "check that the gene symbols use the same nomenclature."); return(NULL)
  }
  message("published_peak_matrix_from_genes(): ", length(genes), " of ",
          dplyr::n_distinct(peaks$gene), " published peak genes present in this cohort.")
  sub <- thr[, c(id_col, genes), drop = FALSE]
  if (!collapse) return(sub)

  # gene -> peak map; a gene inside two peaks contributes to both (GISTIC peaks can overlap)
  g2p <- peaks |> distinct(peak_id, gene) |> filter(gene %in% genes)
  long <- sub |>
    tidyr::pivot_longer(-all_of(id_col), names_to = "gene", values_to = "v") |>
    inner_join(g2p, by = "gene", relationship = "many-to-many") |>
    group_by(.data[[id_col]], peak_id) |>
    summarise(v = mean(suppressWarnings(as.numeric(v)), na.rm = TRUE), .groups = "drop")
  covered <- dplyr::n_distinct(long$peak_id)
  if (covered < dplyr::n_distinct(peaks$peak_id))
    message("  ", covered, " of ", dplyr::n_distinct(peaks$peak_id),
            " published peaks have >=1 gene in this cohort; the rest are filled 0.")
  long |>
    tidyr::complete(!!rlang::sym(id_col), peak_id = unique(peaks$peak_id), fill = list(v = 0)) |>
    pivot_wider(names_from = peak_id, values_from = v, values_fill = 0)
}

# --- GISTIC2 input prep: per-sample .seg -> one combined seg file ------------
# GISTIC2 wants a SINGLE segmentation file with columns
#   Sample  Chromosome  Start  End  Num_Markers  Seg.CN(log2)
# — exactly the DRAGEN .seg content, just concatenated across samples with the
# Sample column set. This reads every data/seg/*.seg (auto-detecting columns, incl.
# the marker-count column), writes the combined file to attend_cnv$seg$gistic_seg_out,
# and returns its path. Run GISTIC on it OUTSIDE R via code/run_gistic.sh; its outputs
# go to data/gistic/ where find_gistic_files() / report 06 picks them up. seg.mean is
# already log2, so NO transform is applied (set value_is_log2=FALSE upstream if not).
# ids: optional character vector of sample ids — when non-NULL, only those samples
# (id-stripped the same way as the combined table) are written, enabling per-group
# GISTIC runs (report 06); returns NULL if ids matches no sample.
# Returns NULL (with a message) when data/seg/ is absent/empty — knit-safe.
write_gistic_seg <- function(cnv = attend_cnv, out = NULL, ids = NULL) {
  sc     <- cnv$seg
  folder <- here("data", sc$dir)
  out    <- out %||% here(sc$gistic_seg_out)
  if (!dir_exists(folder)) { message("write_gistic_seg(): no folder ", folder); return(NULL) }
  files <- dir_ls(folder, recurse = TRUE, type = "file", glob = "*.seg")
  if (length(files) == 0) { message("write_gistic_seg(): no *.seg under ", folder); return(NULL) }

  pick <- function(nms, cands) { h <- intersect(cands, nms); if (length(h)) h[[1]] else NA_character_ }
  combined <- files |>
    map(function(p) {
      id <- str_remove(path_ext_remove(path_file(p)), sc$id_strip)
      dt <- tryCatch(fread(p, header = TRUE, na.strings = c(".", "", "NA"), fill = TRUE) |> as_tibble(),
                     error = function(e) NULL)
      if (is.null(dt) || nrow(dt) == 0) return(NULL)
      cc <- pick(names(dt), sc$chrom_col); ss <- pick(names(dt), sc$start_col)
      ee <- pick(names(dt), sc$end_col);   vv <- pick(names(dt), sc$value_col)
      nn <- pick(names(dt), sc$num_col)
      if (any(is.na(c(cc, ss, ee, vv)))) {
        message("write_gistic_seg(): unresolved columns in ", path_file(p), " — skipping"); return(NULL)
      }
      tibble(Sample      = id,
             Chromosome  = as.character(dt[[cc]]),
             Start       = suppressWarnings(as.numeric(dt[[ss]])),
             End         = suppressWarnings(as.numeric(dt[[ee]])),
             Num_Markers = if (!is.na(nn)) suppressWarnings(as.integer(dt[[nn]])) else NA_integer_,
             Seg.CN      = suppressWarnings(as.numeric(dt[[vv]])))
    }) |>
    purrr::compact() |>
    list_rbind()

  if (nrow(combined) == 0) { message("write_gistic_seg(): nothing to write"); return(NULL) }
  # GISTIC needs Num_Markers; if a .seg lacked it, approximate by segment length in kb
  # (a monotone proxy — GISTIC only uses it for weighting, not position).
  if (any(is.na(combined$Num_Markers)))
    combined$Num_Markers <- ifelse(is.na(combined$Num_Markers),
                                   pmax(1L, as.integer((combined$End - combined$Start) / 1000)),
                                   combined$Num_Markers)
  # Per-group GISTIC runs (report 06) pass ids=. An ids= vector matching NO sample
  # would otherwise write an empty .seg and GISTIC would run on zero samples.
  # combined$Sample is already id-stripped (built above), so normalise ids the same way.
  if (!is.null(ids)) {
    ids <- sub(cnv$seg$id_strip, "", as.character(ids))
    keep <- combined$Sample %in% ids
    if (!any(keep)) {
      message("write_gistic_seg(): ids= matched no sample — nothing written.")
      return(NULL)
    }
    missing <- setdiff(ids, unique(combined$Sample))
    if (length(missing)) {
      message("write_gistic_seg(): ", length(missing),
              " requested id(s) absent from data/seg/: ",
              paste(utils::head(missing, 5), collapse = ", "))
    }
    combined <- combined[keep, , drop = FALSE]
  }
  dir_create(path_dir(out))
  readr::write_tsv(combined, out)
  message("write_gistic_seg(): wrote ", nrow(combined), " segments from ",
          dplyr::n_distinct(combined$Sample), " samples -> ", out)
  out
}

# --- CNV: ASCETS arm-level outputs (report 06, PREFERRED when present) -------
# ASCETS (Spurr et al., Genome Med 2021) converts segmented copy number (e.g.
# DRAGEN .seg) into arm-level SCNAs + a per-sample aneuploidy score. We read the
# `*_arm_weighted_average_segmeans.txt` matrix (the clustering features) and the
# `*_aneuploidy_scores.txt` table. Drop the ASCETS output files into data/ascets/.
# Returns NULL when absent so report 06 falls back to the ClassifyCNV path.
#
# Orientation: ASCETS writes samples as ROWS, arms as COLUMNS. .ascets_orient()
# transposes automatically if a given version emits the matrix the other way.
.ascets_orient <- function(tab) {
  arm_pat <- "^(chr)?([0-9]{1,2}|X|Y)[pq]$"
  if (any(str_detect(names(tab)[-1], arm_pat))) return(tab)   # arms already columns
  long <- tab |> pivot_longer(-1, names_to = "ID", values_to = "v")
  long |> pivot_wider(names_from = 1, values_from = "v")      # arms were rows -> transpose
}

load_ascets_arm_means <- function(cnv = attend_cnv) {
  ac     <- cnv$ascets
  folder <- here("data", ac$dir)
  if (!dir_exists(folder)) { message("load_ascets_arm_means(): no folder ", folder); return(NULL) }
  f <- dir_ls(folder, recurse = TRUE, type = "file", glob = ac$armmeans_glob)
  if (length(f) == 0) { message("load_ascets_arm_means(): no ", ac$armmeans_glob, " under ", folder); return(NULL) }
  tab <- fread(f[[1]], sep = "\t", header = TRUE, na.strings = c("", "NA", "NC")) |>
    as_tibble() |> .ascets_orient()
  names(tab)[1] <- "ID"
  tab |> mutate(ID = str_remove(as.character(ID), ac$id_strip))
}

load_ascets_aneuploidy <- function(cnv = attend_cnv) {
  ac     <- cnv$ascets
  folder <- here("data", ac$dir)
  if (!dir_exists(folder)) return(NULL)
  f <- dir_ls(folder, recurse = TRUE, type = "file", glob = ac$aneuploidy_glob)
  if (length(f) == 0) { message("load_ascets_aneuploidy(): no ", ac$aneuploidy_glob); return(NULL) }
  tab <- fread(f[[1]], sep = "\t", header = TRUE) |> as_tibble()
  names(tab)[1] <- "ID"
  num <- names(tab)[map_lgl(tab, is.numeric)]
  if (length(num) == 0) { message("load_ascets_aneuploidy(): no numeric score column."); return(NULL) }
  tab |> transmute(ID = str_remove(as.character(ID), ac$id_strip),
                   ascets_aneuploidy = suppressWarnings(as.numeric(.data[[num[1]]])))
}

# ASCETS arm-level CALLS (categorical AMP/DEL/NEUTRAL/NC per arm) — the recurrent-CNV
# source for the report-07 oncoplot annotation. Same orientation + id normalisation as
# load_ascets_arm_means(). Returns NULL when absent (knit-safe). Feed to
# recurrent_arm_calls() (attend_classes.R) to rank arms and build the annotation table.
load_ascets_calls <- function(cnv = attend_cnv) {
  ac     <- cnv$ascets
  folder <- here("data", ac$dir)
  if (!dir_exists(folder)) { message("load_ascets_calls(): no folder ", folder); return(NULL) }
  f <- dir_ls(folder, recurse = TRUE, type = "file", glob = ac$calls_glob)
  if (length(f) == 0) { message("load_ascets_calls(): no ", ac$calls_glob, " under ", folder); return(NULL) }
  tab <- fread(f[[1]], sep = "\t", header = TRUE, na.strings = c("", "NA")) |>
    as_tibble() |> .ascets_orient()
  names(tab)[1] <- "ID"
  tab |> mutate(ID = str_remove(as.character(ID), ac$id_strip))
}

# --- GISTIC thresholded matrix in significant peaks = TCGA's EXACT clustering input ----
# Kandoth et al. Suppl. Methods S2: copy-number clustering was run "on thresholded relative
# copy number data in significantly reoccurring amplifications or deletions regions identified
# by GISTIC 2.0". This builds that feature matrix from GISTIC output in data/gistic/:
#   1. read all_thresholded.by_genes.txt (genes x samples, values in {-2,-1,0,1,2});
#   2. restrict to genes inside the SIGNIFICANT peaks (union of amp_genes/del_genes members);
#   3. return an ID x gene tibble (samples as rows) ready for cluster_arm_matrix().
# If the peak-gene lists are absent, falls back to ALL thresholded genes (still thresholded,
# just not peak-restricted) with a message. Returns NULL when GISTIC output is absent — knit-safe.
# The amp_genes/del_genes files list, per significant peak (columns), the member gene symbols
# stacked below a few header rows (cytoband, q-value, residual q, wide peak boundaries).
.gistic_peak_genes <- function(path) {
  if (is.null(path) || !file.exists(path)) return(character(0))
  dt <- tryCatch(fread(path, header = TRUE, sep = "\t", na.strings = c("", "NA")) |> as_tibble(),
                 error = function(e) NULL)
  if (is.null(dt) || ncol(dt) < 2) return(character(0))
  # gene symbols are the cell values below the ~4 annotation rows, across all peak columns
  vals <- unlist(dt[, -1], use.names = FALSE)
  vals <- trimws(as.character(vals))
  hdr  <- c("cytoband", "q value", "q-value", "residual q value", "wide peak boundaries", "genes in wide peak")
  vals <- vals[!is.na(vals) & vals != "" & !tolower(vals) %in% hdr]
  # drop obvious non-gene tokens (coordinates, q-values, bracketed counts)
  vals <- vals[!grepl("^(chr)?[0-9XY]+[:p q0-9.\\-]", vals) & !grepl("^[0-9.eE+-]+$", vals)]
  unique(str_remove(vals, "\\[.*$"))
}

## Core, path-injectable body of load_gistic_thresholded(): reads all_thresholded.by_genes.txt
## from `folder` directly (no here()/data/ resolution), so it is unit-testable against a
## synthetic tempdir. `cnv` still supplies non-path config (id_strip, peak-gene lookup via
## find_gistic_files()). Attaches attr(out, "feature_pos") = data.frame(feature, cytoband) —
## `feature` uses the SAME make.unique() transform as colnames(m) so Task 1's
## feature_positions(colnames(mat), cytoband = ...) can look values up by matrix column name.
load_gistic_thresholded_at <- function(folder, cnv = attend_cnv) {
  g <- cnv$gistic
  if (!dir_exists(folder)) { message("load_gistic_thresholded(): no folder ", folder); return(NULL) }
  ft <- dir_ls(folder, recurse = TRUE, type = "file", glob = g$all_thresholded_glob)
  if (length(ft) == 0) { message("load_gistic_thresholded(): no ", g$all_thresholded_glob, " under ", folder); return(NULL) }

  thr <- tryCatch(fread(ft[[1]], header = TRUE, sep = "\t", na.strings = c("", "NA")) |> as_tibble(),
                  error = function(e) NULL)
  if (is.null(thr) || ncol(thr) < 3) { message("load_gistic_thresholded(): unreadable ", ft[[1]]); return(NULL) }
  # all_thresholded.by_genes: col1 = Gene Symbol, cols 2-3 = Locus ID / Cytoband, rest = samples.
  gene_col <- names(thr)[1]
  meta     <- names(thr)[2:3]
  samp_cols <- setdiff(names(thr), c(gene_col, meta))
  # guard against a trailing-tab phantom column (see load_gistic_lesions_at)
  samp_cols <- samp_cols[!is.na(samp_cols) & nzchar(samp_cols)]

  # restrict to genes in significant peaks (union of amp + del peak members)
  amp <- .gistic_peak_genes(find_gistic_files(cnv)$amp_genes)
  del <- .gistic_peak_genes(find_gistic_files(cnv)$del_genes)
  peak_genes <- union(amp, del)
  genes_all  <- as.character(thr[[gene_col]])
  keep <- if (length(peak_genes) > 0) genes_all %in% peak_genes else rep(TRUE, length(genes_all))
  if (length(peak_genes) == 0)
    message("load_gistic_thresholded(): no peak-gene lists — clustering ALL thresholded genes (still thresholded, not peak-restricted).")
  else
    message("load_gistic_thresholded(): ", sum(keep), " of ", length(genes_all),
            " genes fall in the ", length(peak_genes), " GISTIC significant-peak genes.")

  sub <- thr[keep, , drop = FALSE]
  m <- t(as.matrix(sapply(sub[samp_cols], function(x) suppressWarnings(as.integer(x)))))  # samples x genes
  colnames(m) <- make.unique(as.character(sub[[gene_col]]))
  # meta[2] is the Cytoband column of all_thresholded.by_genes.txt; feature uses the same
  # make.unique() transform as colnames(m) so the two line up by matrix column name.
  fpos <- data.frame(feature  = make.unique(as.character(sub[[gene_col]])),
                      cytoband = as.character(sub[[meta[2]]]),
                      stringsAsFactors = FALSE)
  out <- tibble::as_tibble(m, rownames = "ID") |>
    mutate(ID = str_remove(ID, cnv$seg$id_strip))
  attr(out, "feature_pos") <- fpos
  out
}

load_gistic_thresholded <- function(cnv = attend_cnv) {
  load_gistic_thresholded_at(here("data", cnv$gistic$dir), cnv)
}

# --- GISTIC genome-wide CONTINUOUS matrix = TCGA's Fig-1a DISPLAY, not its clustering input --
# THE DISTINCTION THIS EXISTS FOR. Kandoth et al. use TWO different matrices, and conflating
# them is what made report 09's heatmap not the paper's heatmap:
#
#   CLUSTERING (Suppl. Methods S2) — "thresholded relative copy number data in significantly
#     reoccurring amplifications or deletions regions identified by GISTIC 2.0", i.e.
#     peak-restricted and discretised to {-2..+2}. That is load_gistic_thresholded(), and it
#     stays exactly as it is: it is the faithful method and GISTIC is required for it.
#
#   FIGURE (Fig. 1a legend) — "SCNAs in each tumour (horizontal axis) plotted by chromosomal
#     location (vertical axis)", i.e. the WHOLE GENOME on a CONTINUOUS colour scale.
#
# Plotting the clustering matrix gives a heatmap whose rows are a few hundred peak genes and
# whose colour scale has five steps: genome-ordered, but not the genomic landscape, and not
# continuous. This loader supplies the figure's matrix from the SAME GISTIC run — no .seg
# re-processing and no second tool — so the report can cluster on peaks and display the
# landscape, which is what the paper does.
#
# all_data_by_genes.txt has the same shape as all_thresholded.by_genes.txt (col 1 = Gene
# Symbol, cols 2-3 = Gene ID / Cytoband, rest = samples) but carries CONTINUOUS copy number
# and is NOT peak-restricted. Values are kept as doubles, never coerced to integer, and no
# peak filter is applied — both would undo the point.
#
# Returns an ID x gene tibble with attr "feature_pos" = data.frame(feature, cytoband), the
# same contract load_gistic_thresholded_at() has, so feature_positions(colnames(mat),
# cytoband = ...) resolves rows to chromosome + band ordinal identically. NULL when the file
# is absent, so a report falls back to the clustering matrix and still knits.
load_gistic_continuous_at <- function(folder, cnv = attend_cnv) {
  g <- cnv$gistic
  if (!dir_exists(folder)) { message("load_gistic_continuous(): no folder ", folder); return(NULL) }
  fd <- dir_ls(folder, recurse = TRUE, type = "file", glob = g$all_data_glob)
  if (length(fd) == 0) {
    message("load_gistic_continuous(): no ", g$all_data_glob, " under ", folder,
            " — the Fig-1a display matrix is unavailable; the caller should fall back to the ",
            "thresholded clustering matrix.")
    return(NULL)
  }
  dat <- tryCatch(fread(fd[[1]], header = TRUE, sep = "\t", na.strings = c("", "NA")) |> as_tibble(),
                  error = function(e) NULL)
  if (is.null(dat) || ncol(dat) < 3) { message("load_gistic_continuous(): unreadable ", fd[[1]]); return(NULL) }

  gene_col  <- names(dat)[1]
  meta      <- names(dat)[2:3]
  samp_cols <- setdiff(names(dat), c(gene_col, meta))
  # guard against a trailing-tab phantom column (see load_gistic_lesions_at)
  samp_cols <- samp_cols[!is.na(samp_cols) & nzchar(samp_cols)]
  if (!length(samp_cols)) { message("load_gistic_continuous(): no sample columns in ", fd[[1]]); return(NULL) }

  # as.numeric, NOT as.integer: these are continuous log2-ratio-like values and integer
  # coercion would silently rebuild the discretised matrix this function exists to avoid.
  m <- t(as.matrix(sapply(dat[samp_cols], function(x) suppressWarnings(as.numeric(x)))))
  colnames(m) <- make.unique(as.character(dat[[gene_col]]))
  fpos <- data.frame(feature  = make.unique(as.character(dat[[gene_col]])),
                     cytoband = as.character(dat[[meta[2]]]),
                     stringsAsFactors = FALSE)
  message("load_gistic_continuous(): ", nrow(m), " samples x ", ncol(m),
          " genes genome-wide (continuous) for the Fig-1a display.")
  out <- tibble::as_tibble(m, rownames = "ID") |>
    mutate(ID = str_remove(ID, cnv$seg$id_strip))
  attr(out, "feature_pos") <- fpos
  out
}

load_gistic_continuous <- function(cnv = attend_cnv) {
  load_gistic_continuous_at(here("data", cnv$gistic$dir), cnv)
}

# Locate GISTIC2 output files under data/gistic/ for the report-07 FOCAL overlay
# (paper-standard gene-level recurrent CNV, TCGA nature12113 style). Returns a named
# list(all_lesions, amp_genes, del_genes, scores) of file paths (each NULL if absent), or
# NULL when the folder itself is missing — report 06 then draws mutation+arm only. The
# scores.gistic file is REQUIRED by maftools::readGistic (read.maf forwards the lesion/gene
# files to it and hard-stops without scores), so report 06 only takes the GISTIC overlay
# branch when all four are present. Knit-safe.
find_gistic_files <- function(cnv = attend_cnv) {
  g      <- cnv$gistic
  folder <- here("data", g$dir)
  if (!dir_exists(folder)) return(NULL)
  pick <- function(glb) {
    f <- dir_ls(folder, recurse = TRUE, type = "file", glob = glb)
    if (length(f) == 0) NULL else as.character(f[[1]])
  }
  list(all_lesions = pick(g$all_lesions_glob),
       amp_genes   = pick(g$amp_genes_glob),
       del_genes   = pick(g$del_genes_glob),
       scores      = pick(g$scores_glob))
}

# --- GISTIC all_lesions -> peak x sample matrix -------------------------------
# Report 33 passes all_lesions to maftools as an opaque path; report 06 needs the
# PEAK-LEVEL calls. Peak-level (not gene-level) is the correct testing unit: one
# amplicon spanning 300 genes would otherwise become 300 correlated tests.

#' Parse a GISTIC "Wide Peak Limits" string into chrom/start/end.
#' Input looks like "chr8:127735434-128753674(probes 100:200)".
.gistic_parse_limits <- function(x) {
  x  <- sub("\\(.*$", "", trimws(as.character(x)))   # drop "(probes ...)"
  ch <- sub("^chr", "", sub(":.*$", "", x))
  se <- sub("^[^:]*:", "", x)
  data.frame(chrom      = ch,
             wide_start = suppressWarnings(as.numeric(sub("-.*$", "", se))),
             wide_end   = suppressWarnings(as.numeric(sub("^[^-]*-", "", se))),
             stringsAsFactors = FALSE)
}

#' Read all_lesions.conf_*.txt from an explicit folder (path-injectable core).
#'
#' GISTIC lists every peak TWICE: thresholded rows (0/1/2) and actual-copy-change
#' rows whose Unique Name ends in " - CN values". Only the thresholded rows are kept
#' — keeping both would double every peak and inflate all frequencies 2x.
#'
#' Returns list(peaks = data.frame, mat = samples x peaks matrix), or NULL if absent.
load_gistic_lesions_at <- function(folder, cnv = attend_cnv) {
  if (is.null(folder) || !dir.exists(folder)) return(NULL)
  f <- Sys.glob(file.path(folder, cnv$gistic$all_lesions_glob))
  if (!length(f)) return(NULL)

  les <- tryCatch(
    read.delim(f[1], sep = "\t", check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL)
  if (is.null(les) || !nrow(les)) return(NULL)

  meta_cols <- c("Unique Name", "Descriptor", "Wide Peak Limits", "Peak Limits",
                 "Region Limits", "q values",
                 "Residual q values after removing segments shared with higher peaks",
                 "Broad or Focal", "Amplitude Threshold")
  if (!"Unique Name" %in% names(les)) return(NULL)

  keep <- !grepl("CN values", les[["Unique Name"]], fixed = TRUE)
  les  <- les[keep, , drop = FALSE]
  if (!nrow(les)) return(NULL)

  lim <- .gistic_parse_limits(les[["Wide Peak Limits"]])
  peaks <- data.frame(
    peak_id    = trimws(les[["Unique Name"]]),
    descriptor = trimws(les[["Descriptor"]]),
    direction  = ifelse(grepl("^Amp", trimws(les[["Unique Name"]])), "amp", "del"),
    chrom      = lim$chrom,
    wide_start = lim$wide_start,
    wide_end   = lim$wide_end,
    q_value    = suppressWarnings(as.numeric(les[["q values"]])),
    stringsAsFactors = FALSE
  )

  samp_cols <- setdiff(names(les), meta_cols)
  # GISTIC lesions files carry a trailing tab -> read.delim names the dangling field
  # "" (or NA), which survives setdiff() but is not a real sample; base-R `[` rejects
  # an NA/"" column selection ("undefined columns selected"). Keep real names only.
  samp_cols <- samp_cols[!is.na(samp_cols) & nzchar(samp_cols)]
  if (!length(samp_cols)) return(NULL)

  mat <- t(as.matrix(les[, samp_cols, drop = FALSE]))
  mode(mat) <- "numeric"
  colnames(mat) <- peaks$peak_id
  rownames(mat) <- sub(cnv$seg$id_strip, "", rownames(mat))

  bad <- apply(mat, 1, function(r) all(is.na(r)))
  if (any(bad)) {
    message("load_gistic_lesions_at(): dropping ", sum(bad),
            " non-numeric column(s) misread as samples: ",
            paste(rownames(mat)[bad], collapse = ", "))
    mat <- mat[!bad, , drop = FALSE]
  }

  list(peaks = peaks, mat = mat)
}

#' data/gistic/ wrapper around load_gistic_lesions_at().
load_gistic_lesions <- function(cnv = attend_cnv) {
  load_gistic_lesions_at(here("data", cnv$gistic$dir), cnv)
}

#' Do two GISTIC peaks refer to the same locus?
#'
#' Peak boundaries shift between runs because the confidence interval on peak
#' LOCATION widens as sample count falls. Exact-coordinate matching would call
#' every run's peaks distinct, so overlap is the membership rule.
peaks_overlap <- function(a, b, min_bp = 1) {
  if (as.character(a$chrom[1])     != as.character(b$chrom[1]))     return(FALSE)
  if (as.character(a$direction[1]) != as.character(b$direction[1])) return(FALSE)
  ov <- min(a$wide_end[1], b$wide_end[1]) - max(a$wide_start[1], b$wide_start[1])
  isTRUE(ov + 1 >= min_bp)
}

#' Union of significant peaks across several GISTIC run folders.
#'
#' `folders` must be a NAMED character vector; names become the `source` column.
#' Overlapping peaks from different runs share a `union_id`, so a locus found in
#' three runs is one union peak, not three.
gistic_peak_union <- function(folders, cnv = attend_cnv, min_bp = 1) {
  if (is.null(folders) || !length(folders)) return(NULL)
  if (is.null(names(folders))) names(folders) <- basename(folders)

  parts <- lapply(seq_along(folders), function(i) {
    x <- load_gistic_lesions_at(folders[[i]], cnv)
    if (is.null(x)) return(NULL)
    p <- x$peaks
    p$source <- names(folders)[i]
    p
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(NULL)

  all_p <- do.call(rbind, parts)
  all_p <- all_p[order(all_p$chrom, all_p$direction, all_p$wide_start), , drop = FALSE]

  # Single pass: extend the current union interval while the next peak overlaps it.
  uid <- character(nrow(all_p))
  k <- 0L; cur_end <- -Inf; cur_key <- ""
  for (i in seq_len(nrow(all_p))) {
    key <- paste(all_p$chrom[i], all_p$direction[i])
    if (key != cur_key || all_p$wide_start[i] > cur_end + 1 - min_bp) {
      k <- k + 1L; cur_key <- key; cur_end <- all_p$wide_end[i]
    } else {
      cur_end <- max(cur_end, all_p$wide_end[i])
    }
    uid[i] <- sprintf("U%03d_%s_%s", k, all_p$chrom[i], all_p$direction[i])
  }
  all_p$union_id <- uid
  rownames(all_p) <- NULL
  all_p
}
