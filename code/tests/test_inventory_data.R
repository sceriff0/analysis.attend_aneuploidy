# Smoke-test the PURE logic of inventory_data.R against synthetic fixtures.
# Evaluates only the helper section (everything before the first banner() call),
# so nothing touches disk and data.table/readxl are not needed.
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(purrr); library(tidyr); library(readr); library(fs); library(tibble); library(here)
})
lines <- readLines("code/inventory_data.R")
end   <- grep('^banner\\("0', lines)[1] - 1
start <- grep("^INV <- list\\(", lines)[1]   # skip the library() block; deps loaded above
eval(parse(text = paste(lines[start:end], collapse = "\n")), envir = globalenv())

ok <- 0; bad <- 0
chk <- function(label, got, want) {
  pass <- isTRUE(all.equal(got, want))
  cat(sprintf("  [%s] %-48s got: %s\n", if (pass) "PASS" else "FAIL", label,
              paste(utils::head(as.character(got), 6), collapse = ",")))
  if (pass) ok <<- ok + 1 else { bad <<- bad + 1; cat("        want: ", paste(as.character(want), collapse=","), "\n") }
}

cat("\n-- norm_id (barcode suffix drift the WES loaders paper over)\n")
chk("strips -1TAD104_tumor_only", norm_id("ATT-001-1TAD104_tumor_only"), "ATT-001")
chk("strips _tumor_only",          norm_id("ATT-002_tumor_only"),         "ATT-002")
chk("strips -1TAD104",             norm_id("ATT-003-1TAD104"),            "ATT-003")
chk("strips file extension",       norm_id("ATT-004.seg"),                "ATT-004")
chk("uppercases + trims",          norm_id("  att-005 "),                 "ATT-005")
chk("leaves clean ids alone",      norm_id("ATT-006"),                    "ATT-006")

cat("\n-- col_key (case/separator-insensitive column matching)\n")
col_key <- function(x) str_replace_all(str_to_lower(x), "[^a-z0-9]+", "")
chk("space vs underscore vs dot collapse",
    unique(col_key(c("Unique Subject Identifier","Unique_Subject_Identifier","unique.subject.identifier"))),
    "uniquesubjectidentifier")

cat("\n-- id_flag (which columns are identifiers)\n")
chk("name token wins",   id_flag("Tumor_Sample_Barcode", c("a","a","b")), "name")
chk("shape only",        id_flag("weird_col", paste0("x", 1:50)),         "shape")
chk("both",              id_flag("image_id", paste0("IMG", 1:50)),        "name+shape")
chk("neither (category)",id_flag("MMR_STATUS", rep(c("Deficient","Intact"), 25)), NA_character_)
chk("long free text not an id",
    id_flag("comment", paste0(strrep("z", 60), 1:50)),                     NA_character_)

cat("\n-- profile_columns\n")
df <- tibble(ID = c("A","B","C"), MMR = c("Deficient","Intact",NA), n = c(1,2,3))
p  <- profile_columns(df, "fixture")
chk("one row per column", nrow(p), 3L)
chk("missing counted",    p$n_missing[p$column == "MMR"], 1L)
chk("distinct counted",   p$n_distinct[p$column == "MMR"], 2L)
chk("id detected",        p$id_kind[p$column == "ID"], "name+shape")
chk("low-cardinality prints full vocabulary",
    p$examples[p$column == "MMR"], "{Deficient | Intact}")

cat("\n-- extract_ids\n")
e <- extract_ids(tibble(ID = c("ATT-1_tumor_only","ATT-1_tumor_only","ATT-2")), "fx", "ID")
chk("dedup on raw value", nrow(e), 2L)
chk("normalised alongside raw", sort(e$id_norm), c("ATT-1","ATT-2"))

cat("\n-- overlap math (the join-discovery core)\n")
ids <- bind_rows(
  tibble(source="gianlu", column="TUMOR_BARCODE", id_raw=c("B1","B2","B3"),        id_norm=norm_id(c("B1","B2","B3"))),
  tibble(source="aneu",   column="sequencing_id", id_raw=c("B1-1TAD104","B2_tumor_only"), id_norm=norm_id(c("B1-1TAD104","B2_tumor_only"))),
  tibble(source="orphan", column="ID",            id_raw=c("Z9"),                  id_norm=norm_id("Z9"))
) |> mutate(key = paste0(source, " :: ", column))
uniq <- distinct(ids, key, id_norm); sets <- split(uniq$id_norm, uniq$key); sets <- sets[lengths(sets) > 0]
nm <- names(sets)
ov <- expand_grid(a = nm, b = nm) |> filter(a < b) |>
  mutate(n_a = unname(lengths(sets[a])), n_b = unname(lengths(sets[b])),
         n_shared = unname(map2_int(sets[a], sets[b], ~ length(intersect(.x, .y)))),
         pct_of_smaller = round(100 * n_shared / pmin(n_a, n_b), 1))
chk("all unordered pairs", nrow(ov), 3L)
chk("suffix drift does NOT hide the join", ov$n_shared[ov$n_shared > 0], 2L)
chk("pct of smaller side", ov$pct_of_smaller[ov$n_shared > 0], 100)
linked <- ov |> filter(n_shared > 0)
chk("isolated source found", setdiff(names(sets), union(linked$a, linked$b)), "orphan :: ID")

cat("\n-- empty-input guards (a missing file must not crash the run)\n")
chk("profile_columns(NULL)", nrow(profile_columns(NULL, "x")), 0L)
chk("extract_ids(NULL)",     nrow(extract_ids(NULL, "x", "ID")), 0L)
empty_overlap <- tibble(a=character(), b=character(), n_a=integer(), n_b=integer(),
                        n_shared=integer(), pct_of_smaller=double())
chk("empty overlap still filterable", nrow(filter(empty_overlap, n_shared > 0)), 0L)

cat("\n-- read_table_any dependency guard\n")
tf <- tempfile(fileext = ".csv"); writeLines(c("ID", "A\"B"), tf)
r <- read_table_any(tf)
if (requireNamespace("data.table", quietly = TRUE)) {
  chk("single-column stray quote still reads", is_err(r), FALSE)
} else {
  chk("missing data.table reports the real cause", as.character(r), "data.table not installed")
}
unlink(tf)

cat("\n-- safe() guard (regression: a warning must not cause re-evaluation)\n")
warn_then_value <- function() { warning("noisy but fine"); 42 }
chk("warning muffled, value returned", safe(warn_then_value()), 42)

eval_count <- 0
counted <- function() { eval_count <<- eval_count + 1; warning("w"); eval_count }
invisible(safe(counted()))
chk("expression evaluated EXACTLY once", eval_count, 1)

chk("error captured as probe_err", is_err(safe(stop("boom"))), TRUE)
chk("error message preserved",     as.character(safe(stop("boom"))), "boom")

# The exact shape that escaped the old guard and dropped the session into the
# debugger: fread warns, the handler re-ran it, the re-run errored outside the
# tryCatch. Warn-then-error must now be captured like any other error.
warn_then_error <- function() { warning("first a warning"); stop("then an error") }
chk("warn-then-error is CAPTURED, not escaped", is_err(safe(warn_then_error())), TRUE)

chk("safe_tbl degrades to zero rows", nrow(safe_tbl(stop("unreadable"))), 0L)
chk("safe_tbl passes tables through", nrow(safe_tbl(tibble(a = 1:3))), 3L)

cat(sprintf("\n==== %d passed, %d failed ====\n", ok, bad))
if (bad) quit(status = 1)
