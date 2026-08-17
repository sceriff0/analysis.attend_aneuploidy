# Every inline `r ...` in a report must reference something that report can actually see.
#
# This exists because of a real failure: lifting justification prose out of report 08 into
# 00-methods.Rmd carried an inline `r k_cnv` with it, and k_cnv is defined in report 08's
# setup chunk. Nothing caught it — the file parsed, the chunk bodies parsed, the style rules
# passed — and it only surfaced as "object 'k_cnv' not found" partway through a cluster
# build, after the master had already been rebuilt.
#
# Static, base R, no evaluation: a symbol counts as visible when it is assigned somewhere in
# the report itself, or assigned at top level in one of the code/*.R files THAT REPORT
# SOURCES. Resolving against all of code/ would be wrong — that is exactly the false
# negative that hid this bug, since code/diagnose_unclassified.R happens to set k_cnv and
# no report sources it.

analysis_dir <- file.path("analysis")
# TWO digits, not "^0[0-9]": the glob used to be ^0[0-9][-_] and a report numbered 10 or
# above would have been SILENTLY UNCHECKED — no contract, no derivation notes, no line cap,
# and this suite would still have printed OK. The hyphen rename already came close to
# blinding this pattern once; widening it is cheaper than noticing the gap downstream.
reports <- sort(list.files(analysis_dir, pattern = "^[0-9]{2}[-_].*\\.Rmd$", full.names = TRUE))
stopifnot(length(reports) > 0)

# Base/utility names that appear in inline expressions and are never assigned.
BASE <- c("nrow", "ncol", "length", "sum", "mean", "round", "signif", "paste", "paste0",
          "format", "if", "else", "is.null", "is.na", "max", "min", "seq", "rev", "c",
          "sprintf", "nlevels", "levels", "unique", "sort", "TRUE", "FALSE", "NULL", "NA",
          "function", "return", "collapse", "sep", "dplyr", "n_distinct", "and", "or")

top_level_assigns <- function(txt) {
  m <- regmatches(txt, gregexpr("(?m)^[A-Za-z_.][A-Za-z0-9_.]*\\s*(<-|=)", txt, perl = TRUE))[[1]]
  unique(trimws(sub("\\s*(<-|=)$", "", m)))
}

# Which code/*.R a report sources, following one level of nesting (build_master.R sources
# the loaders, so a report that sources it can see what they define).
sourced_files <- function(txt) {
  hits <- regmatches(txt, gregexpr('source\\(here\\("code", "[^"]+"\\)\\)', txt))[[1]]
  files <- sub('.*"code", "([^"]+)".*', "\\1", hits)
  seen <- character(0)
  repeat {
    new <- setdiff(files, seen)
    if (!length(new)) break
    seen <- c(seen, new)
    for (f in new) {
      p <- file.path("code", f)
      if (!file.exists(p)) next
      s <- paste(readLines(p, warn = FALSE), collapse = "\n")
      nested <- regmatches(s, gregexpr('source\\(here::here\\("code", "[^"]+"\\)\\)|source\\(here\\("code", "[^"]+"\\)\\)', s))[[1]]
      files <- c(files, sub('.*"code", "([^"]+)".*', "\\1", nested))
    }
    files <- unique(files)
  }
  seen
}

fail <- character(0)

for (f in reports) {
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")

  visible <- c(BASE, top_level_assigns(txt))
  for (cf in sourced_files(txt)) {
    p <- file.path("code", cf)
    if (file.exists(p))
      visible <- c(visible, top_level_assigns(paste(readLines(p, warn = FALSE), collapse = "\n")))
  }
  visible <- unique(visible)

  inl <- regmatches(txt, gregexpr("`r [^`]+`", txt))[[1]]
  for (e in inl) {
    expr <- sub("^`r ", "", sub("`$", "", e))
    # Strip numeric literals first, scientific notation included: without this, 1e6 yields
    # a bogus identifier "e6".
    bare <- gsub("[0-9]+\\.?[0-9]*([eE][+-]?[0-9]+)?", " ", expr)
    # Root symbols only: drop anything after $ or ( so field and argument names are ignored.
    roots <- regmatches(bare, gregexpr("[A-Za-z_.][A-Za-z0-9_.]*", bare))[[1]]
    # a name immediately preceded by $ is a field, not a binding
    fields <- regmatches(bare, gregexpr("\\$\\s*[A-Za-z_.][A-Za-z0-9_.]*", bare))[[1]]
    fields <- trimws(sub("^\\$", "", fields))
    roots <- setdiff(roots, fields)
    miss <- setdiff(roots, visible)
    if (length(miss)) {
      fail <- c(fail, paste0(basename(f), ": inline `r ", expr, "` references ",
                             paste(miss, collapse = ", "),
                             ", which this report neither defines nor sources"))
    }
  }
}

if (length(fail)) {
  cat("test_rmd_inline_r: ", length(fail), " unresolved inline reference(s)\n", sep = "")
  cat(paste0("  - ", fail, collapse = "\n"), "\n", sep = "")
} else {
  cat("test_rmd_inline_r: OK (", length(reports), " reports)\n", sep = "")
}
stopifnot(length(fail) == 0L)
