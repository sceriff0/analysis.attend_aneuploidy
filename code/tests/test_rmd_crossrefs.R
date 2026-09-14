# Cross-reference integrity for analysis/*.Rmd.
#
# Base R and parse-free, like test_rmd_style.R: it reads the reports as text, so it runs on a
# machine with no tidyverse and no cohort data exactly as it does on the cluster.
#
# It exists because the 2026-08-13 renumbering moved the FILES and left the PROSE behind, and
# nothing in the suite reads prose. test_rmd_inline_r.R resolves inline `r ...` against the
# source graph; test_plot_style.R reads code lines only. The narrative layer had no guard at
# all, and by 2026-09-14 it had drifted in five files at once:
#
#   * 05-response.Rmd referred to "report 05" EIGHT times, inside report 05, every one of
#     them meaning report 04 (05 used to be 06).
#   * 09-tcga2013-fig1a.Rmd referred to "report 09" TEN times, inside report 09, every one
#     meaning report 07 — including a setup comment reading "Using the SAME pair as report 09
#     is deliberate".
#   * 00-methods.Rmd pointed three parameters at the wrong reports.
#   * 03-covariates-by-mmr.Rmd carried `[05](04-aneuploidy-mmrd.html)` and `[08](06-tmb.html)`:
#     the link LABELS survived a renumbering that repointed their targets, so the site told
#     the reader to go to 05 and sent them to 04.
#
# Four rules:
#   [1] a report never refers to itself by number — say "this report"
#   [2] every "report NN" names a report that exists
#   [3] a link whose LABEL starts with a number agrees with the number in its target
#   [4] every internal .html link resolves to a real analysis/*.Rmd

reports <- sort(list.files(file.path("analysis"), pattern = "^[0-9]{2}[-_].*\\.Rmd$",
                           full.names = TRUE))
stopifnot(length(reports) > 0)

# The report numbers that exist, as two-character strings ("00".."11").
nums <- substr(basename(reports), 1, 2)
# basename without extension, for link-target resolution
stems <- sub("\\.Rmd$", "", basename(reports))

fail <- character(0)
note <- function(...) fail <<- c(fail, paste0(...))

# All matches of a regex, with the 1-based line number each sits on.
hits <- function(lines, re) {
  out <- list()
  for (i in seq_along(lines)) {
    m <- regmatches(lines[i], gregexpr(re, lines[i], perl = TRUE))[[1]]
    if (length(m)) out[[length(out) + 1L]] <- data.frame(n = i, txt = m, stringsAsFactors = FALSE)
  }
  if (length(out)) do.call(rbind, out) else data.frame(n = integer(0), txt = character(0))
}

for (f in reports) {
  lines <- readLines(f, warn = FALSE)
  self  <- substr(basename(f), 1, 2)

  # ---- [1] and [2]: "report NN" / "reports NN" in prose OR in a comment ----------------
  h <- hits(lines, "\\b[Rr]eports?\\s+([0-9]{2})\\b")
  if (nrow(h)) {
    ref <- sub("^.*?([0-9]{2})$", "\\1", h$txt)
    bad <- ref == self
    if (any(bad))
      note(basename(f), ": refers to ITSELF as \"report ", self, "\" at line(s) ",
           paste(unique(h$n[bad]), collapse = ", "),
           " — write \"this report\". A self-reference by number is the signature of a ",
           "renumbering that moved the file and left the prose behind.")
    unknown <- !(ref %in% nums) & !bad
    if (any(unknown))
      note(basename(f), ": names report(s) ", paste(unique(ref[unknown]), collapse = ", "),
           " at line(s) ", paste(unique(h$n[unknown]), collapse = ", "),
           " — no such report exists. Reports are ", nums[1], "-", nums[length(nums)], ".")
  }

  # ---- [3] and [4]: internal links ------------------------------------------------------
  lh <- hits(lines, "\\[[^]]*\\]\\([0-9]{2}-[a-z0-9-]+\\.html[^)]*\\)")
  if (nrow(lh)) {
    for (k in seq_len(nrow(lh))) {
      lab <- sub("^\\[([^]]*)\\].*$", "\\1", lh$txt[k])
      tgt <- sub("^.*\\(([0-9]{2}-[a-z0-9-]+)\\.html.*$", "\\1", lh$txt[k])
      if (!tgt %in% stems)
        note(basename(f), ":", lh$n[k], ": link target ", tgt, ".html has no ",
             "analysis/", tgt, ".Rmd.")
      # A label that STARTS with digits is naming the report number; anything else
      # (a title, a section name) is free text and not checked.
      labn <- regmatches(lab, regexpr("^[0-9]{1,2}", lab))
      if (length(labn) && sprintf("%02d", as.integer(labn)) != substr(tgt, 1, 2))
        note(basename(f), ":", lh$n[k], ": link reads \"", lab, "\" but points at ", tgt,
             " — the label survived a renumbering that repointed the target.")
    }
  }
}

if (length(fail)) {
  cat("test_rmd_crossrefs: FAIL\n"); cat(paste0("  - ", fail, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("test_rmd_crossrefs: ALL PASS (", length(reports), " reports checked)\n", sep = "")
