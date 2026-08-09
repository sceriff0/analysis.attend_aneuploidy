#!/usr/bin/env Rscript
# =============================================================================
# inventory_data.R — WHAT IS ACTUALLY IN data/ ?
# =============================================================================
# Ground truth on every raw input before refactoring the loading layer: every
# source file, every Excel sheet, every column, and every ID space — WITHOUT
# trusting what the loaders claim to read.
#
#     Rscript code/inventory_data.R 2>&1 | tee attend_inventory.txt
#
# Two independent halves, each guarded so one failure never kills the run:
#   RAW    — reads the files on disk directly (loaders NOT involved)
#   LOADER — runs every load_*() and profiles what it emits
# The delta between them is the refactor surface: columns present in the file
# but dropped by the loader, files read by more than one loader, and ID spaces
# that turn out to be the same space under normalisation.
#
# Console output is STRUCTURE ONLY — column names, counts, and a few example
# values. Full ID lists are written to output/inventory/ids.tsv on disk and are
# never printed, so the console log stays safe to paste back.
#
# Machine-readable products (output/inventory/):
#   sources.tsv         one row per source file/sheet: path, format, rows, cols
#   columns.tsv         one row per column: type, missingness, cardinality, examples
#   column_sources.tsv  column name -> which sources carry it (redundancy map)
#   ids.tsv             every unique ID, raw + normalised, per source/column
#   id_overlap.tsv      pairwise intersection of every ID space
#   loader_columns.tsv  what each load_*() emits (LOADER half)
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(fs); library(tidyverse); library(data.table)
})

INV <- list(
  out_dir      = here("output", "inventory"),
  n_preview    = 5000,   # rows read for COLUMN profiling (ID extraction reads all rows)
  dir_sample_n = 3,      # files profiled per directory-source (all files still contribute IDs)
  dir_id_max   = 500,    # cap on files opened for ID extraction inside one directory-source
  max_examples = 6,      # example values shown for a high-cardinality column
  cat_max      = 20,     # <= this many distinct values -> print ALL levels (it is a vocabulary)
  str_trunc    = 60,     # truncate example values to this width
  do_loaders   = TRUE,   # run the LOADER half
  do_ihc       = FALSE   # load_ihc_data() walks the whole FlowPath tree (slow) — opt in
)

dir_create(INV$out_dir, recurse = TRUE)

# --- plumbing ---------------------------------------------------------------
banner <- function(txt) cat("\n", strrep("=", 78), "\n== ", txt, "\n", strrep("=", 78), "\n", sep = "")
# Evaluate ONCE: muffle warnings in place via a calling handler, catch errors via
# tryCatch. (A `warning=` handler in tryCatch cannot re-run the expression — the
# retry would execute outside the tryCatch, so an error there escapes the guard
# and kills the run.)
safe    <- function(expr) {
  tryCatch(
    withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) structure(conditionMessage(e), class = "probe_err")
  )
}
is_err  <- function(x) inherits(x, "probe_err")
# Same guard, but for expressions that must yield a table: failure degrades to
# zero rows so one unreadable source never aborts the inventory.
safe_tbl <- function(expr) { out <- safe(expr); if (is_err(out) || is.null(out)) tibble() else out }
show_df <- function(x, n = Inf) {
  x <- as.data.frame(x)
  if (!nrow(x)) return(invisible(cat("  (none)\n")))
  print(head(x, n), row.names = FALSE)
}
wtsv <- function(x, name) {
  if (is.null(x) || !nrow(x)) return(invisible(NULL))
  readr::write_tsv(x, path(INV$out_dir, paste0(name, ".tsv")))
  cat("  wrote", path(INV$out_dir, paste0(name, ".tsv")), "-", nrow(x), "rows\n")
}

# Delimiter sniff: whichever separator appears most often in the header line.
sniff_delim <- function(path) {
  ln <- tryCatch(readLines(path, n = 50, warn = FALSE), error = function(e) character())
  ln <- ln[!startsWith(ln, "#")]
  if (!length(ln)) return("\t")
  hdr <- ln[[1]]
  counts <- c("\t" = str_count(hdr, fixed("\t")), "," = str_count(hdr, fixed(",")),
              ";"  = str_count(hdr, fixed(";")),  "|" = str_count(hdr, fixed("|")))
  if (max(counts) == 0) return("\t")
  names(counts)[which.max(counts)]
}

# How many leading '#' comment lines before the real header (MAF / cBioPortal style).
header_skip <- function(path) {
  ln <- tryCatch(readLines(path, n = 200, warn = FALSE), error = function(e) character())
  if (!length(ln)) return(0L)
  i <- which(!startsWith(ln, "#"))[1]
  if (is.na(i)) 0L else as.integer(i - 1L)
}

# Generic tabular reader. `select` limits to named columns (the cheap ID pass);
# `n_max = Inf` reads everything.
read_table_any <- function(path, n_max = INV$n_preview, select = NULL) {
  ext <- tolower(path_ext(path))
  if (ext == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) return(structure("arrow not installed", class = "probe_err"))
    return(as_tibble(arrow::read_parquet(path)))
  }
  # Check up front, else a missing data.table would send every file down the
  # unparsed fallback below and hide the real cause behind bogus "UNPARSED" rows.
  if (!requireNamespace("data.table", quietly = TRUE))
    return(structure("data.table not installed", class = "probe_err"))
  # Quoted CSVs need quote handling; genomics flat files routinely carry stray
  # quotes. Try the likely mode first, then fall back to quote-off: fread aborts
  # with "Single column input contains invalid quotes" when a one-column file has
  # an unbalanced quote, and that file is still perfectly readable with quote="".
  q_first <- if (ext %in% c("maf", "seg", "tsv", "txt")) "" else "\""
  q_modes <- unique(c(q_first, ""))
  base_args <- list(input = path, sep = sniff_delim(path), header = TRUE, skip = header_skip(path),
                    fill = TRUE, na.strings = c("", "NA", "."), showProgress = FALSE,
                    data.table = FALSE)
  if (is.finite(n_max)) base_args$nrows <- n_max
  if (!is.null(select)) base_args$select <- select

  out <- NULL
  for (qm in q_modes) {
    out <- safe(do.call(data.table::fread, c(base_args, list(quote = qm))))
    if (!is_err(out)) return(as_tibble(out))
  }
  # Last resort: treat it as a single unparsed column so the source is still
  # inventoried (row count + the header text) instead of vanishing from the report.
  ln <- safe(readLines(path, warn = FALSE))
  if (is_err(ln) || !length(ln)) return(out)   # `out` carries the fread error message
  as_tibble(data.frame(unparsed_line = ln[-1], stringsAsFactors = FALSE)) |>
    (\(d) { names(d) <- paste0("UNPARSED[", substr(ln[1], 1, 40), "]"); d })()
}

# --- ID heuristics ----------------------------------------------------------
# A column is an ID if its NAME carries an identifier token, or its SHAPE looks
# like a key (near-unique short strings). Shape-only hits are kept but flagged,
# because that rule also catches free-text columns.
ID_TOKENS <- c("id", "ids", "pid", "barcode", "barcodes", "sample", "samples", "subject",
               "patient", "case", "seq", "image", "identifier", "name", "tumor")

id_flag <- function(nm, x) {
  tok     <- str_split(str_to_lower(nm), "[^a-z0-9]+")[[1]]
  by_name <- any(tok %in% ID_TOKENS)
  ch      <- as.character(x); ch <- ch[!is.na(ch)]
  n <- length(ch); u <- length(unique(ch))
  by_shape <- n > 0 && u > 1 && u / n > 0.9 && mean(nchar(ch)) < 40
  if (by_name && by_shape) "name+shape" else if (by_name) "name" else if (by_shape) "shape" else NA_character_
}

# The suffix drift the WES loaders already paper over (-1TAD104, _tumor_only).
# Normalising here is what reveals that two "different" ID spaces are one space.
norm_id <- function(x) {
  x <- as.character(x)
  x <- str_remove(x, "-1TAD104_tumor_only$")
  x <- str_remove(x, "_tumor_only$")
  x <- str_remove(x, "-1TAD104$")
  x <- str_remove(x, "\\.(maf|seg|tsv|csv|txt)$")
  str_to_upper(str_trim(x))
}

trunc_str <- function(x) ifelse(nchar(x) > INV$str_trunc, paste0(substr(x, 1, INV$str_trunc), "~"), x)

fmt_examples <- function(x) {
  u <- unique(as.character(x[!is.na(x)]))
  if (!length(u)) return("(all missing)")
  if (length(u) <= INV$cat_max) paste0("{", paste(trunc_str(u), collapse = " | "), "}")
  else paste0(paste(trunc_str(head(u, INV$max_examples)), collapse = " | "), " ...")
}

profile_columns <- function(df, source, sheet = NA_character_) {
  if (is.null(df) || is_err(df) || !ncol(df)) return(tibble())
  tibble(
    source     = source,
    sheet      = sheet,
    position   = seq_along(df),
    column     = names(df),
    type       = unname(map_chr(df, ~ class(.x)[1])),
    n_rows     = nrow(df),
    n_missing  = unname(map_int(df, ~ as.integer(sum(is.na(.x) | (!is.na(.x) & trimws(as.character(.x)) == ""))))),
    n_distinct = unname(map_int(df, ~ as.integer(length(unique(as.character(.x[!is.na(.x)])))))),
    id_kind    = unname(map2_chr(names(df), df, id_flag)),
    examples   = unname(map_chr(df, fmt_examples))
  )
}

extract_ids <- function(df, source, cols, sheet = NA_character_) {
  if (is.null(df) || is_err(df) || !length(cols)) return(tibble())
  cols <- intersect(cols, names(df))
  map_dfr(cols, function(cc) {
    v <- unique(as.character(df[[cc]])); v <- v[!is.na(v) & nzchar(trimws(v))]
    if (!length(v)) return(tibble())
    tibble(source = source, sheet = sheet, column = cc, id_raw = v, id_norm = norm_id(v))
  })
}

# =============================================================================
banner("0. Environment")
cat("R version:", R.version.string, "\n")
cat("project root (here):", here(), "\n")
cat("data/ exists:", dir_exists(here("data")), " | output/clean_data exists:",
    dir_exists(here("output", "clean_data")), "\n")
for (p in c("arrow", "readxl", "data.table", "tidyverse"))
  cat(sprintf("  [%s] %s\n", ifelse(requireNamespace(p, quietly = TRUE), "x", " "), p))

# =============================================================================
banner("1. Source discovery   [Q: what files are actually on disk?]")

# Files the loaders name explicitly. `kind`:
#   table = tabular; excel = multi-sheet workbook; lines = one record per line, not tabular
file_specs <- tribble(
  ~source,            ~path,                                       ~kind,
  "attend_xlsx",      here("data", "attend_data.xlsx"),            "excel",
  "gianlu_clinical",  here("data", "clinical_gianlu_data.csv"),    "table",
  "attend_barcodes",  here("data", "attend_barcodes.csv"),         "table",
  "aneu_scores",      here("data", "aneu_scores.csv"),             "table",
  "hrd_scores",       here("data", "hrd_scores.csv"),              "table",
  "msi_scores",       here("data", "msi_scores.csv"),              "lines",
  "tmb_scores",       here("data", "tmb_scores.csv"),              "lines"
)

dir_specs <- tribble(
  ~source,                ~dir,                                ~glob,
  "variant_annotations",  here("data", "variant_annotations"), "*",
  "seg",                  here("data", "seg"),                 "*.seg",
  "cnv_annotations",      here("data", "cnv_annotations"),     "*.cnv.annotated.tsv",
  "ascets",               here("data", "ascets"),              "*",
  "gistic",               here("data", "gistic"),              "*",
  "tcga_2013",            here("data", "tcga_2013"),           "*",
  "tcga",                 here("data", "tcga"),                "*",
  "nassar",               here("data", "nassar"),              "*",
  "clean_data",           here("output", "clean_data"),        "*"
)

file_specs <- file_specs |> mutate(exists = file_exists(path))
dir_specs  <- dir_specs  |> mutate(exists = dir_exists(dir))

cat("\n-- named files:\n");  show_df(file_specs |> transmute(source, kind, exists, path = as.character(path)))
cat("\n-- directories:\n");  show_df(dir_specs  |> transmute(source, glob, exists, dir = as.character(dir)))

# Anything sitting in data/ that no spec above claims — the unknown unknowns.
known <- c(as.character(file_specs$path))
stray <- if (dir_exists(here("data"))) {
  setdiff(as.character(dir_ls(here("data"), type = "file")), known)
} else character()
cat("\n-- unclaimed files directly in data/ (not referenced by any loader spec):\n")
if (length(stray)) cat(paste0("     ", path_file(stray), collapse = "\n"), "\n") else cat("     (none)\n")

# =============================================================================
banner("2. RAW profile — columns and IDs, straight from the files")

col_rows <- list(); id_rows <- list(); src_rows <- list()

add_source <- function(source, path, format, n_rows, n_cols, note = NA_character_) {
  src_rows[[length(src_rows) + 1]] <<- tibble(
    source = source, path = as.character(path), format = format,
    n_rows = n_rows, n_cols = n_cols, note = note)
}

profile_one <- function(df, source, path, format, sheet = NA_character_, note = NA_character_) {
  if (is_err(df)) { cat(sprintf("  %-24s READ ERROR: %s\n", source, df))
                    add_source(source, path, format, NA, NA, paste("READ ERROR:", df)); return(invisible()) }
  prof <- profile_column_block(df, source, sheet)
  cat(sprintf("  %-24s %6d x %3d\n", source, nrow(df), ncol(df)))
  add_source(source, path, format, nrow(df), ncol(df), note)
  invisible(prof)
}

# Column profiling + a full-file second pass for the ID columns only.
profile_column_block <- function(df, source, sheet = NA_character_) {
  prof <- safe_tbl(profile_columns(df, source, sheet))
  col_rows[[length(col_rows) + 1]] <<- prof
  prof
}

# --- 2a. named tabular files ------------------------------------------------
cat("\n-- tabular files (preview read for shape, full read for IDs):\n")
for (i in which(file_specs$kind == "table" & file_specs$exists)) {
  s <- file_specs$source[i]; p <- file_specs$path[i]
  prev <- read_table_any(p, n_max = INV$n_preview)
  prof <- profile_one(prev, s, p, path_ext(p))
  if (is.null(prof) || !length(prof) || is_err(prev)) next
  idc <- prof$column[!is.na(prof$id_kind)]
  full <- if (length(idc)) read_table_any(p, n_max = Inf, select = idc) else NULL
  if (!is.null(full) && !is_err(full))
    id_rows[[length(id_rows) + 1]] <- safe_tbl(extract_ids(full, s, idc))
}

# --- 2b. the Excel workbook, every sheet ------------------------------------
cat("\n-- Excel workbook sheets:\n")
xi <- which(file_specs$kind == "excel" & file_specs$exists)
if (!length(xi)) cat("  attend_data.xlsx not present\n") else
if (!requireNamespace("readxl", quietly = TRUE)) cat("  readxl not installed\n") else {
  xp     <- file_specs$path[xi[1]]
  sheets <- safe(readxl::excel_sheets(xp))
  if (is_err(sheets)) cat("  ERROR listing sheets:", sheets, "\n") else {
    cat("  sheets:", paste(sheets, collapse = " | "), "\n\n")
    for (sh in sheets) {
      df <- safe(readxl::read_excel(xp, sheet = sh, guess_max = 10000, .name_repair = "unique"))
      src <- paste0("attend_xlsx:", sh)
      prof <- profile_one(df, src, xp, "xlsx", sheet = sh)
      if (is_err(df) || is.null(prof)) next
      idc <- prof$column[!is.na(prof$id_kind)]
      id_rows[[length(id_rows) + 1]] <- safe_tbl(extract_ids(df, src, idc, sheet = sh))
    }
  }
}

# --- 2c. line-oriented files (msi / tmb) ------------------------------------
cat("\n-- line-oriented files (NOT tabular — shown as raw structure):\n")
for (i in which(file_specs$kind == "lines" & file_specs$exists)) {
  s <- file_specs$source[i]; p <- file_specs$path[i]
  ln <- safe(readLines(p, warn = FALSE))
  if (is_err(ln)) { cat(sprintf("  %-24s READ ERROR: %s\n", s, ln)); next }
  cat(sprintf("  %-24s %d lines, %d matching '_tumor_only'\n", s, length(ln),
              sum(str_detect(ln, fixed("_tumor_only")))))
  cat("      first line: ", trunc_str(paste0(substr(ln[1], 1, 200))), "\n", sep = "")
  add_source(s, p, "lines", length(ln), NA, "line-oriented; parsed by regex in load_wes_results.R")
  # The barcode is the leading token before the _tumor_only marker.
  hits <- str_match(ln, "^(.+?)_+tumor_only")[, 2]
  hits <- unique(hits[!is.na(hits)])
  if (length(hits))
    id_rows[[length(id_rows) + 1]] <-
      tibble(source = s, sheet = NA_character_, column = "(regex: leading token)",
             id_raw = hits, id_norm = norm_id(hits))
}

# --- 2d. directory sources --------------------------------------------------
# Profile a SAMPLE of files for shape (they share a schema), but harvest IDs from
# ALL of them — plus the filename stems, which for per-sample files ARE the ID space.
cat("\n-- directory sources:\n")
for (i in which(dir_specs$exists)) {
  s <- dir_specs$source[i]; d <- dir_specs$dir[i]; g <- dir_specs$glob[i]
  files <- safe(dir_ls(d, recurse = TRUE, type = "file", glob = g))
  if (is_err(files) || !length(files)) { cat(sprintf("  %-24s (empty)\n", s)); next }
  by_ext <- split(files, tolower(path_ext(files)))
  cat(sprintf("  %-24s %d files: %s\n", s, length(files),
              paste(sprintf("%s x%d", names(by_ext), lengths(by_ext)), collapse = ", ")))

  # filename stems as an ID space
  stems <- unique(path_ext_remove(path_file(files)))
  id_rows[[length(id_rows) + 1]] <-
    tibble(source = paste0(s, ":filenames"), sheet = NA_character_, column = "(filename stem)",
           id_raw = stems, id_norm = norm_id(stems))
  add_source(paste0(s, ":filenames"), d, "dir", length(files), NA, "IDs derived from file names")

  for (ext in names(by_ext)) {
    fs_ext <- by_ext[[ext]]
    samp   <- head(fs_ext, INV$dir_sample_n)
    sub    <- paste0(s, ":*.", ext)
    prev   <- read_table_any(samp[[1]], n_max = INV$n_preview)
    prof   <- profile_one(prev, sub, samp[[1]], ext,
                          note = sprintf("schema sampled from %s (%d files of this type)",
                                         path_file(samp[[1]]), length(fs_ext)))
    if (is_err(prev) || is.null(prof)) next
    # schema drift check across the sampled files
    if (length(samp) > 1) {
      sigs <- map_chr(samp, function(f) {
        p <- read_table_any(f, n_max = 50)
        if (is_err(p)) "ERROR" else paste(names(p), collapse = ",")
      })
      if (length(unique(sigs)) > 1)
        cat("      !! schema DRIFT across sampled files of this type — columns differ\n")
    }
    idc <- prof$column[!is.na(prof$id_kind)]
    if (!length(idc)) next
    for (f in head(fs_ext, INV$dir_id_max)) {
      one <- read_table_any(f, n_max = Inf, select = idc)
      if (!is_err(one)) id_rows[[length(id_rows) + 1]] <- safe_tbl(extract_ids(one, sub, idc))
    }
    if (length(fs_ext) > INV$dir_id_max)
      cat(sprintf("      (ID harvest capped at %d of %d files — raise INV$dir_id_max)\n",
                  INV$dir_id_max, length(fs_ext)))
  }
}

columns <- bind_rows(col_rows)
ids     <- bind_rows(id_rows) |> distinct(source, column, id_raw, .keep_all = TRUE)
sources <- bind_rows(src_rows)

# =============================================================================
banner("3. Column listing per source   [Q: what does each file actually carry?]")
for (s in unique(columns$source)) {
  blk <- columns |> filter(source == s)
  cat("\n--", s, sprintf("(%d columns)\n", nrow(blk)))
  show_df(blk |> transmute(position, column, type, n_missing, n_distinct,
                           id = coalesce(id_kind, ""), examples))
}

# =============================================================================
banner("4. Column redundancy   [Q: which columns exist in more than one source?]")
# Names are matched case- and separator-insensitively: `Unique Subject Identifier`,
# `Unique_Subject_Identifier` and `unique.subject.identifier` are one column.
col_key <- function(x) str_replace_all(str_to_lower(x), "[^a-z0-9]+", "")
col_src <- columns |>
  mutate(key = col_key(column)) |>
  group_by(key) |>
  summarise(n_sources = n_distinct(source),
            spellings = paste(sort(unique(column)), collapse = " | "),
            sources   = paste(sort(unique(source)), collapse = " | "),
            .groups = "drop") |>
  arrange(desc(n_sources), key)

cat("\n-- shared columns (present in >= 2 sources) — candidates for de-duplication:\n")
show_df(col_src |> filter(n_sources >= 2) |> transmute(n_sources, spellings, sources))

cat("\n-- columns unique to exactly one source — the real payload of that file:\n")
show_df(col_src |> filter(n_sources == 1) |> transmute(spellings, sources), n = 200)

# =============================================================================
banner("5. ID spaces   [Q: how many distinct ID spaces are there, really?]")
id_summary <- ids |>
  group_by(source, column) |>
  summarise(n_unique_raw = n_distinct(id_raw), n_unique_norm = n_distinct(id_norm),
            example = paste(head(sort(unique(id_raw)), 3), collapse = " | "), .groups = "drop") |>
  arrange(desc(n_unique_norm))
show_df(id_summary)

banner("6. ID overlap matrix   [Q: which sources can be joined, and on what?]")
# Overlap is computed on NORMALISED ids, so barcode-suffix drift does not hide a join.
ids  <- ids |> mutate(key = paste0(source, " :: ", column))
uniq <- distinct(ids, key, id_norm)
sets <- split(uniq$id_norm, uniq$key)
sets <- sets[lengths(sets) > 0]

# Fixed shape even when nothing overlaps, so the filters below cannot error.
empty_overlap <- tibble(a = character(), b = character(), n_a = integer(),
                        n_b = integer(), n_shared = integer(), pct_of_smaller = double())
overlap <- if (length(sets) < 2) empty_overlap else {
  nm <- names(sets)
  tidyr::expand_grid(a = nm, b = nm) |>
    filter(a < b) |>
    mutate(n_a      = unname(lengths(sets[a])),
           n_b      = unname(lengths(sets[b])),
           n_shared = unname(map2_int(sets[a], sets[b], ~ length(intersect(.x, .y)))),
           pct_of_smaller = round(100 * n_shared / pmin(n_a, n_b), 1)) |>
    arrange(desc(n_shared))
}
linked <- overlap |> filter(n_shared > 0)
cat("\n-- pairs that share at least one ID (sorted by overlap):\n")
show_df(linked, n = 60)
cat("\n-- ID spaces with NO overlap with anything else (isolated \u2014 cannot be joined):\n")
isolated <- setdiff(names(sets), union(linked$a, linked$b))
if (length(isolated)) cat(paste0("     ", isolated, collapse = "\n"), "\n") else cat("     (none)\n")

# =============================================================================
banner("7. LOADER half   [Q: what does each loader keep, and what does it drop?]")
if (!INV$do_loaders) cat("  (disabled: INV$do_loaders = FALSE)\n") else {
  src <- function(f) tryCatch({ source(here("code", paste0(f, ".R"))); TRUE },
                              error = function(e) { cat("!! could not source", f, ":", conditionMessage(e), "\n"); FALSE })
  invisible(lapply(c("attend_classes", "attend_harmonise", "load_clinical",
                     "load_phenotypes", "load_wes_results"), src))

  loaders <- list(
    clinical_data        = "load_clinical_data",
    imaging_data         = "load_imaging_data",
    gianlu_clinical_data = "load_gianlu_clinical_data",
    aneuploidy_scores    = "load_aneuploidy_scores",
    hrd_scores           = "load_hrd_scores",
    msi_scores           = "load_msi_scores",
    tmb_scores           = "load_tmb_scores",
    maf_data             = "load_maf_data",
    maf_tmb              = "load_maf_tmb",
    cnv_data             = "load_cnv_data",
    seg_data             = "load_seg_data",
    ascets_arm_means     = "load_ascets_arm_means",
    ascets_aneuploidy    = "load_ascets_aneuploidy"
  )
  if (INV$do_ihc) loaders$ihc_data <- "load_ihc_data"

  loader_rows <- list()
  for (nm in names(loaders)) {
    fn <- loaders[[nm]]
    if (!exists(fn, mode = "function")) { cat(sprintf("  %-22s NOT DEFINED\n", nm)); next }
    out <- safe(get(fn)())
    if (is_err(out)) { cat(sprintf("  %-22s ERROR: %s\n", nm, out)); next }
    if (is.null(out) || !is.data.frame(out)) { cat(sprintf("  %-22s returned %s\n", nm, class(out)[1])); next }
    cat(sprintf("  %-22s %6d x %3d : %s\n", nm, nrow(out), ncol(out), paste(names(out), collapse = ", ")))
    loader_rows[[nm]] <- safe_tbl(profile_columns(out, paste0("loader:", nm)))
  }
  loader_cols <- bind_rows(loader_rows)

  cat("\n-- DROP REPORT: columns in the raw source that the loader does not emit\n")
  # Only pairs where we profiled both halves can be compared.
  pairs <- tribble(
    ~loader,                 ~raw,
    "clinical_data",         "attend_xlsx:ClinData",
    "imaging_data",          "attend_xlsx:Imaging (2)",
    "gianlu_clinical_data",  "gianlu_clinical",
    "aneuploidy_scores",     "aneu_scores",
    "hrd_scores",            "hrd_scores"
  )
  for (k in seq_len(nrow(pairs))) {
    lc <- loader_cols |> filter(source == paste0("loader:", pairs$loader[k]))
    rc <- columns     |> filter(source == pairs$raw[k])
    if (!nrow(lc) || !nrow(rc)) next
    dropped <- setdiff(col_key(rc$column), col_key(lc$column))
    added   <- setdiff(col_key(lc$column), col_key(rc$column))
    cat(sprintf("\n  %s  <-  %s\n", pairs$loader[k], pairs$raw[k]))
    cat(sprintf("    raw %d cols, loader emits %d\n", nrow(rc), nrow(lc)))
    cat("    DROPPED:", if (length(dropped)) paste(rc$column[col_key(rc$column) %in% dropped], collapse = " | ") else "(none)", "\n")
    cat("    RENAMED/DERIVED (in loader, not raw):",
        if (length(added)) paste(lc$column[col_key(lc$column) %in% added], collapse = " | ") else "(none)", "\n")
  }
  wtsv(loader_cols, "loader_columns")
}

# =============================================================================
banner("8. Writing machine-readable inventory")
wtsv(sources, "sources")
wtsv(columns, "columns")
wtsv(col_src, "column_sources")
wtsv(ids |> select(source, sheet, column, id_raw, id_norm), "ids")
wtsv(overlap, "id_overlap")

banner("DONE")
cat("Console log above is structure only and safe to paste back.\n")
cat("Full ID lists are on disk at", as.character(INV$out_dir), "— send ids.tsv / id_overlap.tsv\n")
cat("if you want the join analysis done on real values.\n")
