# Every R chunk in every analysis/*.Rmd must parse.
#
# knitr is not installed on the development machine (CLAUDE.md "Local dev caveat"), so this
# extracts chunk bodies as text and hands them to base parse(). That catches the failure
# mode a reorganisation actually introduces — a chunk cut in half, or a brace left behind —
# without needing the cohort data or any package.
#
# It does NOT evaluate anything, so it says nothing about whether a chunk *works*.
#
# It also checks that chunk LABELS are unique within a file. knitr keys its cache and figure
# filenames on the label, so a duplicate is not a warning — it aborts the knit with
# "Duplicate chunk label". That is invisible to the body-parsing loop above (the label lives
# on the fence line, which is never parsed) and to test_rmd_style.R, so a renumbering or a
# copy-pasted guard chunk reached a cluster build twice before this rule existed.

rmd <- sort(c(list.files("analysis", pattern = "\\.Rmd$", full.names = TRUE),
              list.files(file.path("analysis", "_explainers"), pattern = "\\.Rmd$",
                         full.names = TRUE)))
stopifnot(length(rmd) > 0)

fail <- character(0)

for (f in rmd) {
  lines <- readLines(f, warn = FALSE)
  open <- grep("^```\\{[rR][ ,}]", lines)
  close <- grep("^```\\s*$", lines)
  for (o in open) {
    e <- close[close > o]
    if (!length(e)) { fail <- c(fail, paste0(basename(f), ":", o, ": unterminated chunk")); next }
    e <- e[1]
    if (e - o < 2L) next                       # empty chunk (a child= include)
    body <- lines[(o + 1L):(e - 1L)]
    # Inline `r ...` in prose is not our business; only chunk bodies are parsed.
    res <- tryCatch({ parse(text = paste(body, collapse = "\n")); NULL },
                    error = function(err) conditionMessage(err))
    if (!is.null(res)) {
      fail <- c(fail, paste0(basename(f), ":", o, ": ", sub("\n.*", "", res)))
    }
  }
}

# ---- chunk labels are unique within a file ---------------------------------

# The label is the first comma-separated field of the fence header, and only when it is a
# bare name: ```{r, echo=FALSE} is unlabelled (knitr auto-names it), so it can never collide.
chunk_label <- function(header) {
  inner <- sub("^```\\{[rR][ ,]?", "", sub("\\}\\s*$", "", header))
  first <- trimws(strsplit(inner, ",")[[1]][1])
  if (length(first) == 0 || is.na(first) || grepl("=", first)) "" else first
}

for (f in rmd) {
  lines  <- readLines(f, warn = FALSE)
  open   <- grep("^```\\{[rR][ ,}]", lines)
  labels <- vapply(lines[open], chunk_label, character(1), USE.NAMES = FALSE)
  named  <- nzchar(labels)
  dups   <- unique(labels[named][duplicated(labels[named])])
  for (d in dups) {
    at <- open[named][labels[named] == d]
    fail <- c(fail, paste0(basename(f), ": duplicate chunk label '", d, "' at line(s) ",
                           paste(at, collapse = ", "), " — knitr aborts the build on this"))
  }
}

if (length(fail)) {
  cat("test_rmd_parse: ", length(fail), " violation(s)\n", sep = "")
  cat(paste0("  - ", fail, collapse = "\n"), "\n", sep = "")
} else {
  cat("test_rmd_parse: OK (", length(rmd), " files)\n", sep = "")
}
stopifnot(length(fail) == 0L)
