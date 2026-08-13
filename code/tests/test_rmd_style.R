# Prose-standard enforcement for analysis/*.Rmd (spec §7,
# specs/2026-08-11-report-reorganization-design.md).
#
# Parse-only by design: it reads the .Rmd as text and never evaluates a chunk, so it runs
# on a machine with no tidyverse and no cohort data (see CLAUDE.md "Local dev caveat")
# exactly as it does on the cluster. Base R only.
#
# Six rules:
#   [1] every analysis/0N_*.Rmd opens with a scope contract in its first 20 lines
#   [2] every chunk with visible output has a prose line directly above it
#   [3] no analysis/*.Rmd exceeds MAX_LINES
#   [4] no blockquote outside the scope contract (00-methods excepted)
#   [5] no takeaway/interpretation/notes sections anywhere
#   [6] the "two TMB definitions" table lives in exactly one file, an explainer

# The cap exists to stop a report quietly accumulating unrelated analyses. After the merge
# to nine reports, several are DELIBERATELY multi-part ("Part 1 ... Part 2"), so a tight
# cap would fight the chosen structure rather than protect it. It is kept only as a
# runaway guard; the scope contract is what actually keeps a report honest about its scope.
MAX_LINES <- 800L

# Chunks that exist in every report for mechanical reasons and carry no result of their
# own, so rule [2] would only force a sentence restating the label.
BOILERPLATE <- c("setup", "session-info")

analysis_dir <- file.path("analysis")
reports <- sort(list.files(analysis_dir, pattern = "^0[0-9][-_].*\\.Rmd$", full.names = TRUE))
stopifnot(length(reports) > 0)

# Rules [1]-[4] apply to every numbered report. The scope contract covers 00-methods too;
# the blockquote ban does not, because every justification blockquote is *supposed* to end
# up there.
numbered <- reports[!grepl("^00[-_]", basename(reports))]

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
  n <- length(readLines(f, warn = FALSE))
  if (n > MAX_LINES) note(basename(f), ": ", n, " lines, over the ", MAX_LINES, "-line guard")
}

# ---- [4] no blockquotes outside the scope contract -------------------------

for (f in numbered) {
  lines <- readLines(f, warn = FALSE)
  fenced <- in_any_fence(lines)
  quoted <- which(grepl("^>", lines) & !fenced)
  stray <- setdiff(quoted, contract_span(lines))
  if (length(stray)) {
    note(basename(f), ": blockquote outside the scope contract at line(s) ",
         paste(stray, collapse = ", "), " — justification belongs in 00-methods.Rmd")
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

# The explainer fragments were folded into 00-methods.Rmd when the reports merged: with
# nine reports the same block was rendering on up to eight pages, so reports now LINK to
# 00-methods instead of inlining. The table must therefore appear exactly once, there.
all_rmd <- list.files(analysis_dir, pattern = "\\.Rmd$", full.names = TRUE, recursive = TRUE)
carriers <- Filter(function(f) {
  lines <- readLines(f, warn = FALSE)
  any(grepl("^\\|.*tmb__TMB_SCORE.*\\|", lines))
}, all_rmd)
if (length(carriers) != 1L) {
  note("the two-TMB-definitions table appears in ", length(carriers), " file(s) (",
       paste(basename(carriers), collapse = ", "), "); it belongs in 00-methods.Rmd only")
} else if (basename(carriers) != "00-methods.Rmd") {
  note("the two-TMB-definitions table lives in ", basename(carriers), ", not 00-methods.Rmd")
}

# ---- [7] no report re-inlines a shared explainer ---------------------------

for (f in numbered) {
  lines <- readLines(f, warn = FALSE)
  hits <- grep("child\\s*=\\s*'_explainers/", lines)
  if (length(hits)) {
    note(basename(f), ": inlines an explainer at line(s) ", paste(hits, collapse = ", "),
         " — link to 00-methods instead, or the same text renders on every page")
  }
}

# ---- report ----------------------------------------------------------------

if (length(fail)) {
  cat("test_rmd_style: ", length(fail), " violation(s)\n", sep = "")
  cat(paste0("  - ", fail, collapse = "\n"), "\n", sep = "")
} else {
  cat("test_rmd_style: OK\n")
}
stopifnot(length(fail) == 0L)
