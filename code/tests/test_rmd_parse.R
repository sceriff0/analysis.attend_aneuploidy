# Every R chunk in every analysis/*.Rmd must parse.
#
# knitr is not installed on the development machine (CLAUDE.md "Local dev caveat"), so this
# extracts chunk bodies as text and hands them to base parse(). That catches the failure
# mode a reorganisation actually introduces — a chunk cut in half, or a brace left behind —
# without needing the cohort data or any package.
#
# It does NOT evaluate anything, so it says nothing about whether a chunk *works*.

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

if (length(fail)) {
  cat("test_rmd_parse: ", length(fail), " chunk(s) failed to parse\n", sep = "")
  cat(paste0("  - ", fail, collapse = "\n"), "\n", sep = "")
} else {
  cat("test_rmd_parse: OK (", length(rmd), " files)\n", sep = "")
}
stopifnot(length(fail) == 0L)
