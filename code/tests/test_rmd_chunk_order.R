# Forward references between chunks — a variable read in one chunk and assigned only in a
# LATER one.
#
# WHY THIS EXISTS. The tissue-content panels were placed above the immune section but read
# `imm_metrics`, which that section built further down. Every local check passed: the file
# parses, every chunk parses, the style rules hold, and test_rmd_inline_r.R resolves inline
# `r ...` rather than chunk bodies. It knits in an INTERACTIVE session, because the variable
# is already in the global env from a previous run — so the failure is invisible until a
# clean knit on the cluster, where it dies with "object 'imm_metrics' not found" 59% in.
# That is the same shape as the k_cnv bug: visibility that looks fine in the editor and is
# wrong in a fresh process.
#
# It also reports a SECOND shape, added after `ncol` silently resolved to base::ncol.
# A function body used `ncol` as a value while its own signature never declared it, so R
# looked outward, found the base FUNCTION, and the knit died inside is.finite() with
# "default method not implemented for type 'closure'" -- a message that names neither the
# variable nor the chunk. Nothing else can see this: the name resolves, so it is not a
# missing object, and it is not a forward reference either. What marks it is the USE --
# a base function appearing where a value is expected, never called.
#
# Deliberately conservative — it only reports a symbol that
#   (a) is assigned at the top level of some chunk, AND
#   (b) is read in an EARLIER chunk, AND
#   (c) is not assigned anywhere at or before that read, AND
#   (d) is not visible from the source() graph or the base/attached packages, AND
#   (e) is not a data-mask name — a dplyr/ggplot column, which is a symbol in the AST but
#       never a global. Without (e) the check reports `pid` and `p` on reports 01 and 07,
#       which are column names, and a check that cries wolf gets deleted.
# Anything it cannot resolve is passed, so it under-reports rather than blocking a build on
# a guess. Base R only.

rmd_files <- sort(list.files("analysis", pattern = "^[0-9]{2}[-_].*\\.Rmd$", full.names = TRUE))
stopifnot(length(rmd_files) > 0)

# Symbols any report may use without assigning: everything reachable from the sourced
# scripts, plus the packages the reports attach.
src_syms <- unique(unlist(lapply(
  list.files("code", pattern = "\\.R$", full.names = TRUE), function(f) {
    e <- tryCatch(parse(f), error = function(e) NULL)
    if (is.null(e)) return(character(0))
    out <- character(0)
    for (x in as.list(e))
      if (is.call(x) && length(x) >= 3 &&
          as.character(x[[1]]) %in% c("<-", "=", "<<-") && is.name(x[[2]]))
        out <- c(out, as.character(x[[2]]))
    out
  })))

pkg_syms <- unique(unlist(lapply(
  c("base", "stats", "utils", "graphics", "grDevices", "methods"),
  function(p) tryCatch(ls(getNamespace(p)), error = function(e) character(0)))))

# Data-masking verbs: a bare symbol in one of these is a COLUMN, not a global. Collected
# across the whole file and excluded, rather than resolved per call site — a column name is
# a column name everywhere in one report, and the coarse rule keeps this readable.
MASK_VERBS <- c("after_stat", "stat", "select", "mutate", "transmute", "filter", "arrange", "group_by", "ungroup",
                "summarise", "summarize", "count", "add_count", "distinct", "rename", "pull",
                "aes", "across", "all_of", "any_of", "one_of", "tibble", "tribble", "crossing",
                "case_when", "if_else", "n_distinct", "fct_reorder", "fct_rev", "slice_max",
                "slice_min", "semi_join", "anti_join", "left_join", "inner_join", "full_join",
                "pivot_longer", "pivot_wider", "replace_na", "drop_na", "reframe", "rowwise",
                "starts_with", "ends_with", "contains", "matches", "vars", "desc")

# Symbols used as a VALUE: a bare argument, or one side of an operator. Never a call head.
# Passing a function BY NAME is ordinary R -- vapply(rows, is.null, ...), tapply(x, g, sum).
# Arguments to these are skipped, or every apply-family call reports a false positive.
HOF <- c("where","sapply","lapply","vapply","apply","tapply","mapply","Map","Reduce","Filter",
         "Find","Position","do.call","by","aggregate","outer","rapply","eapply","purrr::map",
         "map","map_dbl","map_chr","map_int","map_lgl","map_dfr","imap","pmap","walk",
         "match.fun","Negate","Vectorize","optim","integrate","uniroot","ave","split")
collect_value_uses <- function(x, out, formals_seen = character(0)) {
  if (!is.call(x)) return(invisible())
  h <- x[[1]]; hn <- if (is.name(h)) as.character(h) else ""
  if (hn == "function") {
    fm <- x[[2]]
    if (!is.null(names(fm))) formals_seen <- c(formals_seen, names(fm))
    collect_value_uses(x[[3]], out, formals_seen); return(invisible())
  }
  if (hn %in% c("::", ":::", "$", "@")) return(invisible())
  # A higher-order call: its arguments may legitimately BE functions. Recurse for nested
  # calls but record no bare names, so `sum` in tapply(x, g, sum) is not reported.
  if (hn %in% HOF) {
    for (a in as.list(x)[-1]) if (is.call(a)) collect_value_uses(a, out, formals_seen)
    return(invisible())
  }
  # unname(): the comparison below is against an UNNAMED one-element list, so a named empty
  # argument would slip past it and reach the binding that errors.
  args <- unname(as.list(x)[-1])
  for (k in seq_along(args)) {
    # An empty argument (`x[, 1]`) is R's missing symbol, and BINDING it raises
    # "argument is missing". Single-bracket subsetting compares it without forcing.
    if (identical(args[k], list(quote(expr = )))) next
    a <- args[[k]]
    if (is.name(a)) {
      n <- as.character(a)
      if (nzchar(n) && !(n %in% formals_seen)) out$add(n)
    } else if (is.call(a)) collect_value_uses(a, out, formals_seen)
  }
  invisible()
}

# Every bare name appearing ANYWHERE inside a masking verb, including inside an expression:
# transmute(gpos = offset + end) makes `end` a column, and collect_masked() only sees
# top-level bare arguments. Used ONLY by the base-shadow check. The forward-reference check
# keeps the narrower set, because this one also swallows real globals referenced inside a
# dplyr verb (filter(pid %in% mmrd_pids)) and would blind it to exactly the bug it exists
# for -- the two checks want different amounts of caution, so they get different sets.
collect_masked_deep <- function(x, out, inside = FALSE) {
  if (is.name(x)) { if (inside) out$add(as.character(x)); return(invisible()) }
  if (!is.call(x)) return(invisible())
  h <- x[[1]]; hn <- if (is.name(h)) as.character(h) else ""
  if (hn %in% c("::", ":::", "$", "@")) return(invisible())
  here <- inside || hn %in% MASK_VERBS
  nms <- names(x); if (here && !is.null(nms)) out$add(nms[nzchar(nms)])
  args <- unname(as.list(x)[-1])
  for (k in seq_along(args)) {
    if (identical(args[k], list(quote(expr = )))) next
    collect_masked_deep(args[[k]], out, here)
  }
  invisible()
}

collect_masked <- function(x, out) {
  if (!is.call(x)) return(invisible())
  h <- x[[1]]; hn <- if (is.name(h)) as.character(h) else ""
  if (hn %in% MASK_VERBS) {
    nms <- names(x); if (!is.null(nms)) out$add(nms[nzchar(nms)])
    for (k in seq_along(x)[-1]) if (is.name(x[[k]])) out$add(as.character(x[[k]]))
  }
  for (k in seq_along(x)) if (is.call(x[[k]])) collect_masked(x[[k]], out)
  invisible()
}

# Recursively collect names READ in an expression, and names ASSIGNED at any depth.
walk <- function(x, reads, writes, local) {
  if (is.name(x)) { n <- as.character(x); if (nzchar(n)) reads$add(n); return(invisible()) }
  if (!is.call(x)) return(invisible())
  h <- x[[1]]
  hn <- if (is.name(h)) as.character(h) else ""
  if (hn %in% c("<-", "=", "<<-") && length(x) >= 3) {
    if (is.name(x[[2]])) writes$add(as.character(x[[2]]))
    else walk(x[[2]], reads, writes, local)
    walk(x[[3]], reads, writes, local); return(invisible())
  }
  if (hn == "function") {                      # formals are local, so are their names
    fm <- x[[2]]
    if (!is.null(names(fm))) for (nm in names(fm)) local$add(nm)
    for (d in as.list(fm)) if (!missing(d) && !identical(d, quote(expr = ))) walk(d, reads, writes, local)
    walk(x[[3]], reads, writes, local); return(invisible())
  }
  if (hn %in% c("$", "@") && length(x) >= 3) { walk(x[[2]], reads, writes, local); return(invisible()) }
  # Arguments are values; the head is a call. Recorded apart so a base function used as a
  # value can be told from one that is called.
  if (hn %in% c("::", ":::")) return(invisible())
  if (hn == "for" && length(x) >= 4) { if (is.name(x[[2]])) local$add(as.character(x[[2]]))
    walk(x[[3]], reads, writes, local); walk(x[[4]], reads, writes, local); return(invisible()) }
  for (k in seq_along(x)) if (k > 1 || !is.name(h)) walk(x[[k]], reads, writes, local)
  if (is.name(h)) reads$add(hn)
  invisible()
}
bag <- function() { v <- character(0)
  list(add = function(n) v <<- c(v, n), get = function() unique(v)) }

fail <- character(0)
for (f in rmd_files) {
  L <- readLines(f, warn = FALSE)
  starts <- grep("^```\\{r", L); if (!length(starts)) next
  fences <- grep("^```\\s*$", L)
  chunks <- list()
  for (st in starts) {
    en <- fences[fences > st][1]; if (is.na(en)) next
    lbl <- sub("^```\\{r[, ]*([^,}]*).*$", "\\1", L[st])
    chunks[[length(chunks) + 1L]] <- list(label = trimws(lbl), line = st,
                                          body = paste(L[(st + 1):(en - 1)], collapse = "\n"))
  }
  # One pass for the file's data-mask vocabulary, before any chunk is judged.
  masked <- bag(); masked_deep <- bag()
  for (ch in chunks) {
    e <- tryCatch(parse(text = ch$body), error = function(e) NULL)
    if (!is.null(e)) for (x in as.list(e)) {
      collect_masked(x, masked); collect_masked_deep(x, masked_deep)
    }
  }
  mask_syms <- masked$get()
  mask_deep <- unique(c(mask_syms, masked_deep$get()))

  assigned <- character(0)                      # everything written at or before this chunk
  per <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    e <- tryCatch(parse(text = chunks[[i]]$body), error = function(e) NULL)
    if (is.null(e)) { per[[i]] <- list(reads = character(0), writes = character(0)); next }
    r <- bag(); w <- bag(); lo <- bag(); vu <- bag()
    for (x in as.list(e)) walk(x, r, w, lo)
    for (x in as.list(e)) collect_value_uses(x, vu)
    per[[i]] <- list(reads = setdiff(r$get(), lo$get()), writes = w$get(),
                     value_uses = setdiff(vu$get(), lo$get()))
  }
  all_writes <- unique(unlist(lapply(per, `[[`, "writes")))
  for (i in seq_along(chunks)) {
    unresolved <- setdiff(per[[i]]$reads,
                          c(assigned, per[[i]]$writes, src_syms, pkg_syms, mask_syms))
    later <- intersect(unresolved, all_writes)  # assigned, but only in a LATER chunk
    # A base function used as a VALUE and never assigned anywhere in the report: the name
    # resolves, so nothing errors until the value reaches something that wants a number.
    shadow <- setdiff(intersect(per[[i]]$value_uses, pkg_syms),
                      c(assigned, per[[i]]$writes, all_writes, src_syms, mask_deep))
    shadow <- setdiff(shadow, c("T", "F", "pi", "LETTERS", "letters", "month.name",
                                "month.abb", "TRUE", "FALSE", "NULL", "NA", "Inf", "NaN"))
    if (length(shadow))
      fail <- c(fail, sprintf("%s [%s] line %d: uses %s as a VALUE, but it is only a base function",
                              basename(f), chunks[[i]]$label, chunks[[i]]$line,
                              paste(sort(shadow), collapse = ", ")))
    if (length(later))
      fail <- c(fail, sprintf("%s [%s] line %d: reads %s, assigned only in a later chunk",
                              basename(f), chunks[[i]]$label, chunks[[i]]$line,
                              paste(sort(later), collapse = ", ")))
    assigned <- unique(c(assigned, per[[i]]$writes))
  }
}

if (length(fail)) {
  cat("test_rmd_chunk_order:", length(fail), "forward reference(s)\n")
  cat(paste0("  - ", fail, collapse = "\n"), "\n")
  cat("\nA chunk may only read what an EARLIER chunk assigned. This knits in an interactive\n",
      "session and fails on a clean build. Move the assignment above its first use.\n", sep = "")
}
stopifnot(length(fail) == 0L)
cat("test_rmd_chunk_order: OK (", length(rmd_files), " files)\n", sep = "")
