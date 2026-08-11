# Task 6, fix round 1: load_cnv_data()'s attr(out, "n_profiled") must count DISTINCT
# derived sample ids, not raw file count. cc$id_strip (e.g. "-1TAD104|_tumor_only")
# strips MULTIPLE sequencing-id suffixes down to one bare barcode BY DESIGN — the same
# sample can be resequenced/reannotated under a different suffix, or (with
# dir_ls(recurse = TRUE)) the same sample's file can sit under two subfolders. Either
# way, length(files) would overcount samples profiled, inflating segment_pileup()'s
# denominator and deflating every reported recurrence frequency — the opposite-direction
# twin of the bug this branch exists to fix.
#
# The fix extracted the id derivation into .classifycnv_file_id(p, cc) in
# code/load_wes_results.R, used both inside load_cnv_data()'s map() body and for the
# attribute (via dplyr::n_distinct() over ids derived from `files`), so the two paths
# cannot drift.
#
# This test exercises .classifycnv_file_id() ITSELF, not a reimplementation of it, but
# cannot source() code/load_wes_results.R directly: that file's top-level
# library(here)/library(fs)/library(data.table) are not installed on this machine (see
# the Task 6 report). Instead it parse()s the file (no execution of the library() calls
# or any other top-level code) and eval()s just the one function definition, in an
# environment supplying stringr::str_remove() (installed) and a path_file() stand-in
# (base::basename() — identical to fs::path_file() for a plain path with no trailing
# slash, which is all a file path from dir_ls() ever is).
exprs <- parse(file.path("code", "load_wes_results.R"))
is_named_fn <- function(name) vapply(exprs, function(e) {
  is.call(e) && length(e) >= 2 &&
    identical(e[[1]], as.name("<-")) &&
    identical(e[[2]], as.name(name))
}, logical(1))
is_id_fn <- is_named_fn(".classifycnv_file_id")
stopifnot(sum(is_id_fn) == 1)   # exactly one definition in the file

env <- new.env()
env$str_remove <- stringr::str_remove
env$path_file  <- base::basename
eval(exprs[[which(is_id_fn)]], envir = env)
id_fn <- env$.classifycnv_file_id
stopifnot(is.function(id_fn))

# .classifycnv_n_profiled(files, cc) — the pure helper load_cnv_data() calls for
# attr(out, "n_profiled"), extracted so this behaviour is testable without fs/data.table
# (Task 6 fix round 2: reverting load_wes_results.R's n_distinct() back to length(files)
# used to pass the whole suite silently, because load_cnv_data() cannot run here — this
# test closes that hole by exercising the extracted count itself, not a source-grep).
is_np_fn <- is_named_fn(".classifycnv_n_profiled")
stopifnot(sum(is_np_fn) == 1)
env$n_distinct <- dplyr::n_distinct
eval(exprs[[which(is_np_fn)]], envir = env)
n_profiled_fn <- env$.classifycnv_n_profiled
stopifnot(is.function(n_profiled_fn))

cc <- list(id_strip = "-1TAD104|_tumor_only")

# [1] Two filenames differing only by the id_strip suffix -> the SAME id (the exact
# collision the finding raised: S1-1TAD104.cnv.annotated.tsv and
# S1_tumor_only.cnv.annotated.tsv both derive "S1").
id_a <- id_fn("data/cnv_annotations/S1-1TAD104.cnv.annotated.tsv", cc)
id_b <- id_fn("data/cnv_annotations/S1_tumor_only.cnv.annotated.tsv", cc)
stopifnot(identical(id_a, "S1"))
stopifnot(identical(id_b, "S1"))
stopifnot(identical(id_a, id_b))

# [2] The same collision holds when the two files sit in different subfolders
# (dir_ls(recurse = TRUE) would list both).
id_c <- id_fn("data/cnv_annotations/batch1/S1-1TAD104.cnv.annotated.tsv", cc)
id_d <- id_fn("data/cnv_annotations/batch2/S1_tumor_only.cnv.annotated.tsv", cc)
stopifnot(identical(id_c, "S1"))
stopifnot(identical(id_c, id_d))

# [3] Distinct samples still derive distinct ids (the fix must not over-collapse).
id_s2 <- id_fn("S2-1TAD104.cnv.annotated.tsv", cc)
stopifnot(identical(id_s2, "S2"))
stopifnot(!identical(id_a, id_s2))

# [4] n_distinct() over a 2-file, 1-sample vector gives 1, not 2 -- the exact quantity
# load_cnv_data() now attaches as attr(out, "n_profiled").
ids <- vapply(
  c("S1-1TAD104.cnv.annotated.tsv", "S1_tumor_only.cnv.annotated.tsv"),
  id_fn, character(1), cc = cc
)
stopifnot(length(ids) == 2)
stopifnot(dplyr::n_distinct(ids) == 1)

# [5] .classifycnv_n_profiled(): the SAME 2-file/1-sample collision, exercised through the
# extracted helper itself, not a reimplementation. Reverting load_cnv_data() to
# length(files) would make this return 2, not 1 -- pins the regression the plan flagged.
stopifnot(n_profiled_fn(c("S1-1TAD104.cnv.annotated.tsv",
                          "S1_tumor_only.cnv.annotated.tsv"), cc) == 1)

# [6] Distinct samples still count distinctly through the helper (3 files, 2 samples: S1
# resequenced once, S2 once).
stopifnot(n_profiled_fn(c("S1-1TAD104.cnv.annotated.tsv",
                          "S1_tumor_only.cnv.annotated.tsv",
                          "S2-1TAD104.cnv.annotated.tsv"), cc) == 2)

cat("All .classifycnv_file_id() de-duplication properties verified.\n")
