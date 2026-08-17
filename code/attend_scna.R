# =============================================================================
# attend_scna.R — statistics for report 06 (recurrent SCNAs by aneuploidy x MMR)
#
# Sourced by analysis/07-mutation-and-cna.Rmd (Part 3). Depends on
# attend_classes.R (attend_scna, attend_cnv) being sourced first.
#
# Design constraint: MMRd-high is n=9. Between-group PER-PEAK inference is not
# supported at that size (spec §2.1); the primary endpoint is an AGGREGATED panel
# score. Helpers here are split accordingly — powered aggregates vs explicitly
# exploratory per-peak tables.
# =============================================================================

#' Per-sample segment counts from a folder of .seg files (path-injectable).
seg_counts_at <- function(seg_dir, cnv = attend_cnv) {
  if (is.null(seg_dir) || !dir.exists(seg_dir)) return(NULL)
  f <- list.files(seg_dir, pattern = "\\.seg$", full.names = TRUE, recursive = TRUE)
  if (!length(f)) return(NULL)
  out <- do.call(rbind, lapply(f, function(p) {
    n <- tryCatch(nrow(read.delim(p, sep = "\t", check.names = FALSE)),
                  error = function(e) NA_integer_)
    data.frame(ID = sub(cnv$seg$id_strip, "", tools::file_path_sans_ext(basename(p))),
               n_seg = as.integer(n), stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' Audit which samples GISTIC will silently drop at -maxseg, and whether those
#' exclusions are associated with scna_group.
#'
#' Hypersegmentation correlates with aneuploidy (more real breakpoints -> more
#' segments), so a maxseg cut preferentially removes aneuploidy-high tumours. With
#' MMRd-high at n=9, losing 3 shrinks the group of interest by a third and biases
#' the composite score downward for purely technical reasons.
maxseg_audit <- function(seg_counts, groups, maxseg = attend_scna$maxseg) {
  m <- merge(seg_counts, groups, by = "ID", all.x = TRUE)
  m$excluded <- !is.na(m$n_seg) & m$n_seg > maxseg

  lev <- levels(groups$scna_group)
  tab <- table(scna_group = factor(as.character(m$scna_group), levels = lev),
               excluded   = factor(m$excluded, levels = c(FALSE, TRUE)))

  fp <- if (any(m$excluded) && length(unique(na.omit(as.character(m$scna_group)))) > 1) {
    tryCatch(stats::fisher.test(tab, simulate.p.value = TRUE, B = 20000)$p.value,
             error = function(e) NA_real_)
  } else NA_real_

  list(excluded  = m[m$excluded, c("ID", "n_seg", "scna_group"), drop = FALSE],
       tab       = tab,
       fisher_p  = fp,
       n_before  = nrow(m),
       n_after   = sum(!m$excluded))
}

# --- per-group / leave-one-out .seg writers (inputs for run_gistic.sh MODE=groups|loo) ---
# write_gistic_seg() alone emits ONLY the pooled seg (attend_all_segments.seg); the group
# and LOO runs have no input until these run. All three share ONE barcode->group mapping
# (barcode_scna_group) so a GISTIC group can never disagree with the report's grouping.

#' Barcodes (id-stripped .seg stems), IDENTICAL to seg_counts_at()$ID and to
#' write_gistic_seg()'s per-sample Sample — filenames only, no file contents read.
seg_barcodes <- function(cnv = attend_cnv) {
  seg_dir <- here("data", cnv$seg$dir)
  if (!dir.exists(seg_dir)) return(character(0))
  f <- list.files(seg_dir, pattern = "\\.seg$", full.names = TRUE, recursive = TRUE)
  sub(cnv$seg$id_strip, "", tools::file_path_sans_ext(basename(f)))
}

#' scna_group factor for a vector of BARCODES, via the barcode->pid crosswalk.
#' THE single source of truth report 06's group_of_barcode() delegates to, so the
#' .seg writers and the report cannot assign a barcode to different groups.
barcode_scna_group <- function(barcodes, master, cw,
                               pid_col = attend_cols$pid,
                               lev     = attend_scna$group_levels) {
  if (is.null(master) || is.null(cw) || !nrow(cw))
    return(factor(rep(NA_character_, length(barcodes)), levels = lev))
  pid_v <- cw$pid[match(barcodes, cw$barcode)]
  factor(as.character(master$scna_group[match(pid_v, master[[pid_col]])]), levels = lev)
}

#' .seg filename token for a group level: "MMRd-high" -> "mmrd_high" (matches the
#' attend_<token>.seg / attend_<token>_drop_* names run_gistic.sh globs).
scna_group_token <- function(level) tolower(gsub("-", "_", level))

#' Write one combined attend_<token>.seg per scna_group for MODE=groups. Reuses
#' write_gistic_seg(ids=), so each group's per-file column resolution and Num_Markers
#' handling are byte-for-byte the same path as the pooled seg. Returns written paths.
write_group_segs <- function(master, cw, cnv = attend_cnv, segdir = NULL) {
  if (is.null(segdir)) segdir <- dirname(here(cnv$seg$gistic_seg_out))
  bc  <- seg_barcodes(cnv)
  grp <- barcode_scna_group(bc, master, cw)
  written <- character(0)
  for (level in attend_scna$group_levels) {
    ids <- bc[!is.na(grp) & grp == level]
    out <- file.path(segdir, paste0("attend_", scna_group_token(level), ".seg"))
    if (!length(ids)) { message("write_group_segs(): no samples for ", level, " — skipped."); next }
    res <- write_gistic_seg(cnv = cnv, out = out, ids = ids)
    if (!is.null(res)) written <- c(written, stats::setNames(res, scna_group_token(level)))
  }
  message("write_group_segs(): wrote ", length(written), " group seg(s) -> ", segdir)
  invisible(written)
}

#' Write leave-one-out .seg for one group (default MMRd-high): one
#' attend_<token>_drop_<barcode>.seg per member, each dropping that member, for
#' MODE=loo peak-stability (spec §4.2). Needs >=3 members; returns written paths.
write_loo_segs <- function(master, cw, cnv = attend_cnv, segdir = NULL,
                           group = "MMRd-high") {
  if (is.null(segdir)) segdir <- dirname(here(cnv$seg$gistic_seg_out))
  bc  <- seg_barcodes(cnv)
  grp <- barcode_scna_group(bc, master, cw)
  hi  <- bc[!is.na(grp) & grp == group]
  if (length(hi) < 3) {
    message("write_loo_segs(): only ", length(hi), " ", group,
            " sample(s) — LOO needs >=3; skipped.")
    return(invisible(character(0)))
  }
  token   <- scna_group_token(group)
  written <- character(0)
  for (id in hi) {
    out <- file.path(segdir, paste0("attend_", token, "_drop_", id, ".seg"))
    res <- write_gistic_seg(cnv = cnv, out = out, ids = setdiff(hi, id))
    if (!is.null(res)) written <- c(written, res)
  }
  message("write_loo_segs(): wrote ", length(written), " LOO seg(s) for ",
          group, " (n=", length(hi), ") -> ", segdir)
  invisible(written)
}

#' Map each pre-specified panel locus to a GISTIC peak of the MATCHING direction.
#'
#' Direction is part of the match, not a filter applied afterwards: 8q24 can carry
#' both an amplification and a deletion peak, and only the amplification is the
#' serous-like MYC event.
match_panel_peaks <- function(peaks, panel = attend_scna$panel) {
  hit <- vapply(seq_len(nrow(panel)), function(i) {
    ok <- peaks$direction == panel$direction[i] &
      startsWith(as.character(peaks$descriptor), panel$cytoband[i])
    if (!any(ok)) return(NA_character_)
    # Most significant peak wins when a locus matches several.
    peaks$peak_id[ok][which.min(peaks$q_value[ok])]
  }, character(1))

  data.frame(locus = panel$locus, direction = panel$direction,
             peak_id = hit, stringsAsFactors = FALSE)
}

#' Composite "serous-like" score: fraction of the pre-specified TCGA UCEC panel
#' altered per sample, in each locus's expected direction.
#'
#' Report 15's PRIMARY endpoint. Aggregating is ~4x more powerful than per-peak
#' testing at n=9 (spec §2.1): simulated power 93-100% for moderate effects vs 22%
#' for a single typical peak, uncorrected.
#'
#' Loci with no matching peak are dropped from the denominator rather than scored
#' as unaltered — counting them as zero would deflate every sample equally and
#' shrink the observed effect.
panel_score <- function(mat, peaks, panel = attend_scna$panel) {
  mp  <- match_panel_peaks(peaks, panel)
  mp  <- mp[!is.na(mp$peak_id) & mp$peak_id %in% colnames(mat), , drop = FALSE]
  n   <- nrow(mp)

  if (!n) {
    out <- stats::setNames(rep(NA_real_, nrow(mat)), rownames(mat))
    attr(out, "n_loci") <- 0L
    return(out)
  }

  sub <- mat[, mp$peak_id, drop = FALSE]
  out <- rowSums(sub >= 1, na.rm = TRUE) / n
  names(out) <- rownames(mat)
  attr(out, "n_loci")  <- n
  attr(out, "loci")    <- mp$locus
  out
}

#' Per-peak alteration frequency by group, with WILSON score CIs.
#'
#' [U5] Wilson (prop.test, correct = FALSE) is preferred over Clopper-Pearson at
#' small n — CP is over-conservative and over-wide (Brown, Cai & DasGupta 2001).
#' The interval width is the point at n = 9 (4/9 -> ~17-69%): it is displayed, not
#' hidden. Wilson is base R, so no dependency is added.
.wilson_ci <- function(k, n) {
  if (n <= 0) return(c(NA_real_, NA_real_))
  # prop.test warns on small n / extreme counts; the Wilson interval is still valid.
  suppressWarnings(as.numeric(stats::prop.test(k, n, correct = FALSE)$conf.int))
}

scna_freq_table <- function(mat, group) {
  stopifnot(nrow(mat) == length(group))
  lev <- levels(group)
  if (is.null(lev)) lev <- sort(unique(as.character(group)))

  rows <- list()
  for (g in lev) {
    idx <- which(as.character(group) == g)
    if (!length(idx)) next
    sub <- mat[idx, , drop = FALSE]
    for (p in colnames(mat)) {
      k <- sum(sub[, p] >= 1, na.rm = TRUE); n <- length(idx)
      ci <- .wilson_ci(k, n)
      rows[[length(rows) + 1L]] <- data.frame(
        peak_id = p, scna_group = g, n = n, n_altered = k,
        freq = if (n > 0) k / n else NA_real_,
        ci_lo = ci[1], ci_hi = ci[2], stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  out$scna_group <- factor(out$scna_group, levels = lev)
  rownames(out) <- NULL
  out
}

#' Two-group label-permutation test on a continuous score.
#'
#' Permuting the label preserves the observed burden distribution exactly, which
#' matters here: aneuploidy-high tumours carry more SCNAs by construction, so any
#' test that does not hold the burden fixed is confounded by design.
#'
#' The (1 + k) / (B + 1) correction keeps p strictly positive — an uncorrected 0/B
#' would print as p = 0, a claim n = 9 cannot support.
perm_test_two_group <- function(x, g, B = attend_scna$perm_B, seed = 1) {
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]; g <- droplevels(factor(g[keep]))
  if (nlevels(g) != 2) stop("perm_test_two_group(): need exactly 2 groups, got ", nlevels(g))

  lv   <- levels(g)
  obs  <- mean(x[g == lv[1]]) - mean(x[g == lv[2]])
  set.seed(seed)
  null <- replicate(B, {
    gp <- sample(g)
    mean(x[gp == lv[1]]) - mean(x[gp == lv[2]])
  })

  list(stat = obs,
       p    = (1 + sum(abs(null) >= abs(obs))) / (B + 1),
       B    = B)
}

#' Correct within a testing family. Family A (12 pre-specified loci) uses Holm;
#' Family B (genome-wide discovery) uses BH. They are NEVER pooled — pooling would
#' spend the confirmatory panel's power on discovery multiplicity.
#'
#' `n` defaults to `length(p)` — the FULL length of the input vector, including any
#' NA / untested entries. This is deliberate and differs from `stats::p.adjust()`'s
#' own default: p.adjust drops NA internally before sizing its correction, so calling
#' it directly on a vector with NAs corrects over only the non-NA (tested) subset —
#' the SMALLER, wrong denominator. [U1] requires the LARGER denominator — the full
#' union peak count (Family B) or the full 12-locus panel (Family A) — so that a
#' locus/peak that didn't resolve to a testable p-value still counts toward the
#' family's multiplicity; under-counting would inflate apparent significance. For a
#' p-vector with no NA, `n = length(p)` reproduces the prior (pre-fix) result exactly,
#' so this change is backward-compatible.
family_adjust <- function(p, family = c("A", "B"), cfg = attend_scna, n = length(p)) {
  family <- match.arg(family)
  stats::p.adjust(p, method = if (family == "A") cfg$panel_p_adjust else cfg$fdr_method, n = n)
}

#' Shared / private structure of per-run peak sets, with Jaccard overlap.
#'
#' Layer 3 of the report. At n=9 "is MYC differentially frequent?" is unanswerable,
#' but "do the two groups recurrently alter the same loci?" is — as a set-overlap
#' question rather than a per-peak test.
peak_set_structure <- function(union_tbl) {
  srcs <- unique(union_tbl$source)
  sets <- lapply(srcs, function(s) unique(union_tbl$union_id[union_tbl$source == s]))
  names(sets) <- srcs

  jac <- matrix(NA_real_, length(srcs), length(srcs), dimnames = list(srcs, srcs))
  for (i in seq_along(srcs)) for (j in seq_along(srcs)) {
    a <- sets[[i]]; b <- sets[[j]]
    un <- length(union(a, b))
    jac[i, j] <- if (un == 0) NA_real_ else length(intersect(a, b)) / un
  }

  tally   <- table(unique(union_tbl[, c("source", "union_id")])$union_id)
  private <- unique(union_tbl[union_tbl$union_id %in% names(tally)[tally == 1],
                              c("source", "union_id"), drop = FALSE])
  rownames(private) <- NULL

  list(sets = sets, jaccard = jac, private = private)
}

#' Which reference group does each target sample's CN profile resemble?
#'
#' Layer 2. Stays powered at n=9 because it is a per-sample summary, not a per-peak
#' test: each of the 9 gets a correlation against each reference centroid.
profile_correlation <- function(mat, group, target, refs) {
  gi  <- as.character(group)
  idx <- which(gi == target)
  if (!length(idx)) return(NULL)

  cent <- lapply(refs, function(r) {
    k <- which(gi == r)
    if (!length(k)) return(NULL)
    colMeans(mat[k, , drop = FALSE], na.rm = TRUE)
  })
  names(cent) <- refs
  cent <- cent[!vapply(cent, is.null, logical(1))]
  if (!length(cent)) return(NULL)

  out <- data.frame(ID = rownames(mat)[idx], stringsAsFactors = FALSE)
  for (r in names(cent)) {
    out[[r]] <- vapply(idx, function(i)
      suppressWarnings(stats::cor(mat[i, ], cent[[r]], use = "complete.obs")),
      numeric(1))
  }
  out$best_match <- names(cent)[max.col(as.matrix(out[, names(cent), drop = FALSE]),
                                        ties.method = "first")]
  rownames(out) <- NULL
  out
}

#' Fraction of leave-one-out runs that retain each full-run peak.
#'
#' A peak driven by a single tumour vanishes when that tumour is dropped; one
#' surviving all 9 runs is real recurrence. "Retained" = a significant peak in the
#' LOO run whose wide-limit interval overlaps the full-run peak by >= min_bp, on the
#' same chromosome and in the same direction. Exact-coordinate matching would
#' understate stability, since boundaries shift between runs.
loo_stability <- function(full_peaks, loo_folders, cnv = attend_cnv,
                          min_bp = attend_scna$loo_min_overlap_bp) {
  if (is.null(full_peaks) || !nrow(full_peaks) ||
      is.null(loo_folders) || !length(loo_folders)) return(NULL)

  loo <- lapply(loo_folders, function(f) {
    x <- load_gistic_lesions_at(f, cnv); if (is.null(x)) NULL else x$peaks
  })
  loo <- loo[!vapply(loo, is.null, logical(1))]
  if (!length(loo)) return(NULL)

  n_ret <- vapply(seq_len(nrow(full_peaks)), function(i) {
    sum(vapply(loo, function(lp)
      any(vapply(seq_len(nrow(lp)), function(j)
        peaks_overlap(full_peaks[i, ], lp[j, ], min_bp), logical(1))),
      logical(1)))
  }, integer(1))

  data.frame(peak_id       = full_peaks$peak_id,
             n_loo         = length(loo),
             n_retained    = n_ret,
             retained_frac = n_ret / length(loo),
             stringsAsFactors = FALSE)
}

# --- data-driven panel: selection NESTED inside the permutation --------------
# The Family B counterpart to panel_score(). panel_score() asks "are these tumours
# serous-like in the TCGA sense?" from a panel fixed before the data were seen. This
# asks "is there ANY directional peak set that separates the two groups?" — a
# strictly larger question, answered from the data, and therefore Family B forever.
#
# Correlated with nothing in Family A, and never pooled with it: family_adjust()
# keeps the two families apart precisely so a discovery result cannot spend the
# confirmatory panel's power.

#' Altered-in-the-given-direction, as a logical matrix.
#'
#' `dir` is +1 (amplification) or -1 (deletion) per column. Multiplying by dir
#' flips a deletion's sign, so one `>= 1` test serves both directions: a -2 call at
#' a del-direction locus becomes +2 and counts, while a +2 call there becomes -2 and
#' does not. Same "concordant events only" rule match_panel_peaks() applies to the
#' pre-specified panel — an undirected rule answers "most altered", not "most like
#' the target phenotype".
.directional_altered <- function(mat, idx, dir) {
  sub <- mat[, idx, drop = FALSE]
  sweep(sub, 2, dir, `*`) >= 1
}

#' Choose the n_select most group-discriminating peaks, with their directions.
#'
#' DELIBERATELY label-dependent — that is the whole point, and the reason the result
#' can never be treated as pre-specified. Criterion is the signed difference in mean
#' thresholded call between the two groups: magnitude ranks the peaks, sign fixes
#' each one's direction. One criterion doing both jobs, so there is no second
#' selection step to account for.
select_panel_loci <- function(mat, grp, n_select = attend_scna$select_n) {
  lv <- levels(droplevels(factor(grp)))
  if (length(lv) != 2L) stop("select_panel_loci(): need exactly 2 groups, got ", length(lv))
  d  <- colMeans(mat[grp == lv[1], , drop = FALSE], na.rm = TRUE) -
        colMeans(mat[grp == lv[2], , drop = FALSE], na.rm = TRUE)
  d[!is.finite(d)] <- 0
  # A peak with a zero difference carries no direction, so it is not selectable:
  # sign(0) == 0 would make .directional_altered() test `0 >= 1` for every sample
  # and contribute a constant, silently shrinking the score's range.
  cand <- which(d != 0)
  if (!length(cand)) return(list(idx = integer(0), dir = numeric(0)))
  idx <- cand[order(abs(d[cand]), decreasing = TRUE)][seq_len(min(n_select, length(cand)))]
  list(idx = idx, dir = sign(d[idx]))
}

#' Group-mean difference in data-driven panel score, for one labelling.
.selected_panel_stat <- function(mat, grp, n_select) {
  lv  <- levels(droplevels(factor(grp)))
  sel <- select_panel_loci(mat, grp, n_select)
  if (!length(sel$idx)) return(0)
  s   <- rowMeans(.directional_altered(mat, sel$idx, sel$dir), na.rm = TRUE)
  mean(s[grp == lv[1]], na.rm = TRUE) - mean(s[grp == lv[2]], na.rm = TRUE)
}

#' Permutation test for a data-driven panel, with the selection re-run inside every
#' replicate.
#'
#' ⚠️ THE VALIDITY IS ENTIRELY IN WHERE THE SELECTION HAPPENS. select_panel_loci()
#' is called INSIDE the replicate, on the PERMUTED labels. Hoisting it out — select
#' once, then permute the scores — is the invalid "double dipping" version
#' (Kriegeskorte et al., Nat Neurosci 2009): the panel would be chosen with the real
#' labels, so the permuted scores could never reproduce the optimism the selection
#' introduced, the null would be far too narrow, and p would be near zero on pure
#' noise. That hoist looks like an obvious speed-up and is invisible in a diff.
#' test_nested_selection_permutation.R pins the CALIBRATION (uniform p under the
#' null), which is what actually catches it.
#'
#' ONE-SIDED by construction, not by choice: select_panel_loci() takes each locus's
#' direction from the sign of the group difference, so the selected panel always
#' favours the first group and the statistic cannot be negative. Comparing |null| to
#' |obs|, as perm_test_two_group() does for a fixed panel, would be testing a
#' two-sided hypothesis the statistic cannot express.
perm_test_selected_panel <- function(mat, grp,
                                     n_select = attend_scna$select_n,
                                     B        = attend_scna$select_perm_B,
                                     seed     = 1) {
  keep <- !is.na(grp)
  mat  <- mat[keep, , drop = FALSE]
  grp  <- droplevels(factor(grp[keep]))
  if (nlevels(grp) != 2L)
    stop("perm_test_selected_panel(): need exactly 2 groups, got ", nlevels(grp))

  obs <- .selected_panel_stat(mat, grp, n_select)
  set.seed(seed)
  null <- replicate(B, .selected_panel_stat(mat, sample(grp), n_select))

  sel <- select_panel_loci(mat, grp, n_select)
  # Leave-one-out selection stability. A panel whose membership turns over when one
  # patient is dropped is not a finding, and at n ~ 9 per group that is the likely
  # outcome — so it is reported beside p rather than left for a reader to wonder about.
  loo <- table(unlist(lapply(seq_len(nrow(mat)), function(i)
    colnames(mat)[select_panel_loci(mat[-i, , drop = FALSE], grp[-i], n_select)$idx])))

  list(stat = obs,
       p    = (1 + sum(null >= obs)) / (B + 1),
       B    = B,
       loci = colnames(mat)[sel$idx],
       dir  = ifelse(sel$dir > 0, "amp", "del"),
       loo_frac = as.numeric(loo[colnames(mat)[sel$idx]]) / nrow(mat))
}

#' Global two-group test over the WHOLE peak set, selecting nothing.
#'
#' The complement to both panels: total directional burden per sample (any call of
#' |value| >= 1 at any union peak), through the existing fixed-panel permutation. No
#' selection, so no selection to correct for, and it stays answerable when both panel
#' routes come back null — it asks "is the profile different at all?", which neither
#' a serous-aimed panel nor a discriminating-peak panel can answer.
perm_test_global_burden <- function(mat, grp, B = attend_scna$perm_B, seed = 1) {
  burden <- rowMeans(abs(mat) >= 1, na.rm = TRUE)
  perm_test_two_group(burden, grp, B = B, seed = seed)
}
