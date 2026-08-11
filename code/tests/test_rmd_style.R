# Prose-standard enforcement for analysis/*.Rmd (spec §7,
# specs/2026-08-11-report-reorganization-design.md).
#
# Parse-only by design: it reads the .Rmd as text and never evaluates a chunk, so it runs
# on a machine with no tidyverse and no cohort data (see CLAUDE.md "Local dev caveat")
# exactly as it does on the cluster. Base R only.
#
# Six rules:
#   [1] every analysis/[0-4]*.Rmd opens with a scope contract in its first 20 lines
#   [2] every chunk with visible output has a prose line directly above it
#   [3] no analysis/*.Rmd exceeds MAX_LINES
#   [4] no blockquote in analysis/[1-4]*.Rmd outside the scope contract
#   [5] no takeaway/interpretation/notes sections anywhere
#   [6] the "two TMB definitions" table lives in exactly one file, an explainer

MAX_LINES <- 400L

# Spec section 10 anticipated that the 400-line cap is a first guess and told us to adjust it
# rather than force an unnatural split. Adjusting it globally would licence every report to
# grow, so the exemptions are named one at a time with their own ceiling and their reason.
# A file whose exemption is no longer needed fails too (STALE), so these cannot outlive
# the argument for them.
OVER_CAP <- list(
  # The single destination for every justification in the pipeline (D5). Length is the
  # point: capping it pushes prose back into the reports it was extracted from.
  `00_methods.Rmd` = 500L,
  # Spec section 4 maps 04:358-773 here as one report. It carries four question-groups
  # (composition, immune content, continuous covariates, imaging calls) and is the
  # natural next split if the report set is ever extended past the 15 the spec plans.
  `32_cell_composition.Rmd` = 550L,
  # Named in spec section 10 as the file expected to exceed the cap even after the
  # methods extraction. It does: clustering, k-sweep, enrichment, and the cascade.
  `40_tcga_classification.Rmd` = 650L
)

# Chunks that exist in every report for mechanical reasons and carry no result of their
# own, so rule [2] would only force a sentence restating the label.
BOILERPLATE <- c("setup", "session-info")

analysis_dir <- file.path("analysis")
reports <- sort(list.files(analysis_dir, pattern = "^[0-4][0-9]_.*\\.Rmd$", full.names = TRUE))
stopifnot(length(reports) > 0)

# [1-4] apply to the numbered reports; the scope contract also covers 00_methods.Rmd
# (spec §7 globs [0-4]*), whereas the blockquote ban starts at 1x because every
# justification blockquote is *supposed* to land in 00_methods.
banded <- reports[!grepl("^00_", basename(reports))]

fail <- character(0)
note <- function(...) fail <<- c(fail, paste0(...))

# ---- helpers ---------------------------------------------------------------

# Chunk fences, as (start, end) line pairs. A fence line is ```{r ...} to open and ``` to
# close; anything between them is code, so headers/blockquotes inside chunks (i.e. R
# comments) must not be mistaken for prose by the rules below.
chunk_spans <- function(lines) {
  opens <- grep("^```\\{[rR][ ,}]", lines)
  if (length(opens) == 0) return(list(start = integer(0), end = integer(0)))
  ends <- integer(length(opens))
  for (i in seq_along(opens)) {
    after <- which(lines == "```" & seq_along(lines) > opens[i])
    ends[i] <- if (length(after)) after[1] else length(lines)
  }
  list(start = opens, end = ends)
}

# TRUE for every line that sits inside a fenced block of any language, so prose rules
# only ever look at real prose.
in_any_fence <- function(lines) {
  fence <- grepl("^```", lines)
  inside <- logical(length(lines))
  open <- FALSE
  for (i in seq_along(lines)) {
    if (fence[i]) { inside[i] <- TRUE; open <- !open; next }
    inside[i] <- open
  }
  inside
}

# Prose positions: line numbers that are neither YAML front matter nor inside a fence.
# The scope contract is located against THIS index rather than raw line numbers, because a
# contract that quotes a live `r nrow(dat)` must sit after the setup chunk that builds dat
# (knitr evaluates inline R in document order). "First 20 lines" therefore means the first
# 20 lines a reader actually reads.
prose_index <- function(lines) {
  keep <- !in_any_fence(lines)
  if (length(lines) && trimws(lines[1]) == "---") {
    close_yaml <- which(trimws(lines) == "---")
    close_yaml <- close_yaml[close_yaml > 1L]
    if (length(close_yaml)) keep[1:close_yaml[1]] <- FALSE
  }
  which(keep)
}

# The contract is the contiguous run of blockquote lines, inside the prose head, that
# carries the four required keys.
contract_span <- function(lines) {
  pi <- prose_index(lines)
  head_lines <- utils::head(pi, 20)
  q <- head_lines[grepl("^>", lines[head_lines])]
  if (!length(q)) return(integer(0))
  # extend the run downward: a contract may wrap past the 20th prose line
  last <- max(q)
  while (last + 1L <= length(lines) && grepl("^>", lines[last + 1L])) last <- last + 1L
  min(q):last
}

chunk_label <- function(header) {
  inner <- sub("^```\\{[rR][ ,]?", "", sub("\\}\\s*$", "", header))
  first <- trimws(strsplit(inner, ",")[[1]][1])
  if (length(first) == 0 || is.na(first) || grepl("=", first)) "" else first
}

has_opt_false <- function(header, opt) {
  grepl(paste0(opt, "\\s*=\\s*FALSE"), header)
}

# ---- [1] scope contract ----------------------------------------------------

CONTRACT_KEYS <- c("\\*\\*Question\\.\\*\\*", "\\*\\*Cohort\\.\\*\\*",
                   "\\*\\*Reads\\.\\*\\*", "\\*\\*Out of scope\\.\\*\\*")

for (f in reports) {
  lines <- readLines(f, warn = FALSE)
  span <- contract_span(lines)
  head_txt <- if (length(span)) lines[span] else character(0)
  missing <- CONTRACT_KEYS[!vapply(CONTRACT_KEYS, function(k) any(grepl(k, head_txt)), logical(1))]
  if (length(missing)) {
    note(basename(f), ": scope contract incomplete in the first 20 prose lines (missing ",
         paste(gsub("\\\\|\\*", "", missing), collapse = ", "), ")")
  }
}

# ---- [2] derivation note above every visible-output chunk ------------------

for (f in reports) {
  lines <- readLines(f, warn = FALSE)
  sp <- chunk_spans(lines)
  for (i in seq_along(sp$start)) {
    h <- lines[sp$start[i]]
    lbl <- chunk_label(h)
    # A child= chunk is an inclusion, not a result: the derivation note for the figure it
    # precedes comes after it (spec §5.3 worked example).
    if (grepl("child\\s*=", h)) next
    if (lbl %in% BOILERPLATE) next
    if (has_opt_false(h, "include")) next
    j <- sp$start[i] - 1L
    while (j >= 1L && !nzchar(trimws(lines[j]))) j <- j - 1L
    ok <- j >= 1L &&
      !grepl("^#", lines[j]) &&        # a header is not a derivation note
      !grepl("^```", lines[j]) &&      # nor is the tail of another chunk
      !grepl("^>", lines[j])           # nor the scope contract
    if (!ok) {
      note(basename(f), ":", sp$start[i], ": chunk '", if (nzchar(lbl)) lbl else "<unlabelled>",
           "' has no prose line directly above it")
    }
  }
}

# ---- [3] length cap --------------------------------------------------------

for (f in list.files(analysis_dir, pattern = "\\.Rmd$", full.names = TRUE)) {
  n   <- length(readLines(f, warn = FALSE))
  b   <- basename(f)
  cap <- if (!is.null(OVER_CAP[[b]])) OVER_CAP[[b]] else MAX_LINES
  if (n > cap) {
    note(b, ": ", n, " lines, over its ", cap, "-line cap")
  } else if (!is.null(OVER_CAP[[b]]) && n <= MAX_LINES) {
    note(b, ": ", n, " lines — now under the ", MAX_LINES,
         "-line cap, so its OVER_CAP exemption is STALE and should be deleted")
  }
}

# ---- [4] no blockquotes outside the scope contract -------------------------

for (f in banded) {
  lines <- readLines(f, warn = FALSE)
  fenced <- in_any_fence(lines)
  quoted <- which(grepl("^>", lines) & !fenced)
  stray <- setdiff(quoted, contract_span(lines))
  if (length(stray)) {
    note(basename(f), ": blockquote outside the scope contract at line(s) ",
         paste(stray, collapse = ", "), " — justification belongs in 00_methods.Rmd")
  }
}

# ---- [5] no takeaway sections ----------------------------------------------

BANNED_HEADERS <- "^##+ *(What to take away|Interpretation|Notes) *$"
for (f in list.files(analysis_dir, pattern = "\\.Rmd$", full.names = TRUE)) {
  lines <- readLines(f, warn = FALSE)
  hits <- which(grepl(BANNED_HEADERS, lines) & !in_any_fence(lines))
  if (length(hits)) {
    note(basename(f), ": banned section header at line(s) ", paste(hits, collapse = ", "),
         " — reports show derivations and results (D6)")
  }
}

# ---- [6] the two-TMB-definitions table is written once ---------------------

all_rmd <- list.files(analysis_dir, pattern = "\\.Rmd$", full.names = TRUE, recursive = TRUE)
carriers <- Filter(function(f) {
  lines <- readLines(f, warn = FALSE)
  any(grepl("^\\|.*tmb__TMB_SCORE.*\\|", lines))
}, all_rmd)
if (length(carriers) != 1L) {
  note("the two-TMB-definitions table appears in ", length(carriers), " file(s) (",
       paste(basename(carriers), collapse = ", "), "); it belongs in _explainers/ only")
} else if (basename(dirname(carriers)) != "_explainers") {
  note("the two-TMB-definitions table lives in ", basename(carriers),
       ", not in analysis/_explainers/")
}

# ---- report ----------------------------------------------------------------

if (length(fail)) {
  cat("test_rmd_style: ", length(fail), " violation(s)\n", sep = "")
  cat(paste0("  - ", fail, collapse = "\n"), "\n", sep = "")
} else {
  cat("test_rmd_style: OK\n")
}
stopifnot(length(fail) == 0L)
