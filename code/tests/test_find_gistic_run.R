# find_gistic_files() must resolve to a NAMED run, never to whichever folder sorts first.
#
# WHY THIS EXISTS. data/gistic/ holds one subfolder per GISTIC run — the pooled cohort and
# the four scna_group strata, plus loo/ — and find_gistic_files() used to glob the whole
# tree with recurse = TRUE and take f[[1]]. It resolved to `all/` only because "all" sorts
# before "loo" and "mmr*". Nothing checked that. A run folder named earlier in sort order
# (an "ATTEND_v2/", a "batch1/") would have silently repointed report 07's clustering peaks
# and report 08's focal overlay at a 9-patient stratum — every report still knitting, every
# figure still drawing, every number wrong. That is the failure mode this suite exists for:
# the wrong answer looks exactly like the right one.
#
# Four properties are pinned:
#   [1] with no `dir`, the POOLED run is chosen BY NAME, even when another folder sorts first
#   [2] with `dir`, that exact folder is used and no other run can leak in
#   [3] a flat data/gistic/ with no per-run subfolders still resolves (the legacy layout)
#   [4] a missing folder returns NULL, so a report degrades to a skip rather than erroring
# Base R only, on a synthetic tree — no cohort data, no cluster.

suppressWarnings(suppressMessages(library(fs)))

# Evaluate ONLY find_gistic_files() out of the loader, rather than source()ing it. The file
# opens with library(tidyverse), which this machine does not have (CLAUDE.md "Local dev
# caveat"), and a run-selection check that can only run on the cluster cannot catch a
# run-selection bug before it reaches the cluster — which is the whole point of this file.
# fs and here() are everything this one function touches.
.exprs <- parse(file.path("code", "load_wes_results.R"))
for (.e in as.list(.exprs))
  if (is.call(.e) && length(.e) >= 3 && identical(as.character(.e[[1]]), "<-") &&
      is.name(.e[[2]]) && identical(as.character(.e[[2]]), "find_gistic_files"))
    eval(.e, envir = globalenv())
if (!exists("find_gistic_files", mode = "function"))
  stop("find_gistic_files() not found as a top-level assignment in code/load_wes_results.R")

fail <- character(0)
note <- function(...) fail <<- c(fail, paste0(...))

GFILES <- c("all_lesions.conf_99.txt", "amp_genes.conf_99.txt",
            "del_genes.conf_99.txt", "scores.gistic")
# File CONTENTS are the run name, so every assertion below checks provenance rather than
# merely "a file was found".
mkrun <- function(d, tag) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  for (f in GFILES) writeLines(tag, file.path(d, f))
}
check <- function(got, want, tag) {
  if (is.null(got)) { note("[", tag, "] returned NULL"); return(invisible()) }
  for (k in names(got)) {
    if (is.null(got[[k]])) { note("[", tag, "] ", k, " not found"); next }
    saw <- readLines(got[[k]])[1]
    if (!identical(saw, want))
      note("[", tag, "] ", k, " resolved to run '", saw, "', not '", want, "'")
  }
}

cfg <- list(gistic = list(
  dir = "gistic", pooled_dir = "all",
  all_lesions_glob = "*all_lesions*.txt", amp_genes_glob = "*amp_genes*.txt",
  del_genes_glob   = "*del_genes*.txt",   scores_glob    = "*scores.gistic"))

base <- file.path(tempdir(), paste0("gx_", as.integer(runif(1, 1, 1e6))))
gx   <- file.path(base, "data", "gistic")
# "aaa_first" exists only to sort BEFORE "all" — it is the adversary the old code lost to.
for (r in c("aaa_first", "all", "mmrd_high")) mkrun(file.path(gx, r), r)

# here() is shadowed in globalenv, where find_gistic_files() was evaluated, so the default
# (no-`dir`) branch resolves into the synthetic tree instead of the real repo.
here <- function(...) file.path(base, ...)

# ---- [1] no dir -> the pooled run, by name ---------------------------------
check(find_gistic_files(cfg), "all", "1")

# ---- [2] explicit dir -> that run, and nothing else ------------------------
check(find_gistic_files(cfg, dir = file.path(gx, "mmrd_high")), "mmrd_high", "2")

# ---- [3] a flat tree with no per-run subfolders still works ----------------
flat <- file.path(base, "flat")
mkrun(flat, "flat")
check(find_gistic_files(cfg, dir = flat), "flat", "3")

# ---- [4] a missing folder is NULL, not an error ----------------------------
if (!is.null(find_gistic_files(cfg, dir = file.path(base, "definitely_absent"))))
  note("[4] a missing folder must return NULL so reports degrade to a skip")

unlink(base, recursive = TRUE)

if (length(fail)) {
  cat("test_find_gistic_run: FAIL\n"); cat(paste0("  - ", fail, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("test_find_gistic_run: OK\n")
