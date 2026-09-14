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
#' ⚠️ This does NOT control for SCNA burden. Permuting labels keeps the pooled score
#' distribution but not a burden difference BETWEEN the groups: if MMRp-high is simply
#' more aneuploid than MMRd-high (both sit above one cut), any fixed panel scores
#' higher there for that reason alone. panel_specificity_test() is the burden control.
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
#' Which GISTIC runs exist, and which are USABLE — as a table a report can print.
#'
#' run_gistic.sh guards its group sweep per group (`|| echo WARN ... continuing`), so one
#' group failing leaves the rest intact and a report simply finds nothing where that group
#' should be: a section that renders empty rather than one that says a run is missing.
#' readGistic() additionally needs ALL FOUR of all_lesions / amp_genes / del_genes /
#' scores.gistic and hard-stops without them, so "the folder exists" is not "the run is
#' usable". Both distinctions are reported here rather than inferred from a blank page.
gistic_run_inventory <- function(run_dirs, cnv = attend_cnv) {
  need <- c(cnv$gistic$all_lesions_glob, cnv$gistic$amp_genes_glob,
            cnv$gistic$del_genes_glob,   cnv$gistic$scores_glob)
  do.call(rbind, lapply(names(run_dirs), function(nm) {
    d <- run_dirs[[nm]]
    if (!dir.exists(d))
      return(data.frame(run = nm, files = "0/4", status = "MISSING - no folder",
                        stringsAsFactors = FALSE))
    got <- vapply(need, function(g) length(Sys.glob(file.path(d, g))) > 0, logical(1))
    # A folder holding only gistic_inputs.mat is GISTIC started and died; distinguishing that
    # from "never attempted" is what tells you whether to fix the .seg or the GISTIC call.
    started <- length(list.files(d)) > 0
    data.frame(run = nm, files = paste0(sum(got), "/4"),
               status = if (all(got)) "complete"
                        else if (started) paste0("INCOMPLETE - started, missing ",
                                                 paste(gsub("[*]", "", need[!got]), collapse = ", "))
                        else "EMPTY - folder created, nothing written",
               stringsAsFactors = FALSE)
  }))
}

#' Jaccard is confounded by run size: GISTIC's power grows with n, so a nine-patient
#' run finds few peaks and its Jaccard against a large run is low BY CONSTRUCTION even
#' when every peak it finds is shared. The overlap coefficient |A n B| / min(|A|, |B|)
#' is 1 when the smaller set is a subset of the larger, so it is returned beside Jaccard
#' with each run's peak count, and is the one to read for "does the small run find
#' anything the large one does not".
peak_set_structure <- function(union_tbl) {
  srcs <- unique(union_tbl$source)
  sets <- lapply(srcs, function(s) unique(union_tbl$union_id[union_tbl$source == s]))
  names(sets) <- srcs

  jac <- ovl <- matrix(NA_real_, length(srcs), length(srcs), dimnames = list(srcs, srcs))
  for (i in seq_along(srcs)) for (j in seq_along(srcs)) {
    a <- sets[[i]]; b <- sets[[j]]
    un <- length(union(a, b)); mn <- min(length(a), length(b))
    jac[i, j] <- if (un == 0) NA_real_ else length(intersect(a, b)) / un
    ovl[i, j] <- if (mn == 0) NA_real_ else length(intersect(a, b)) / mn
  }

  tally   <- table(unique(union_tbl[, c("source", "union_id")])$union_id)
  private <- unique(union_tbl[union_tbl$union_id %in% names(tally)[tally == 1],
                              c("source", "union_id"), drop = FALSE])
  rownames(private) <- NULL

  list(sets = sets, jaccard = jac, overlap = ovl,
       n_peaks = vapply(sets, length, integer(1)), private = private)
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
# asks "is there ANY peak set that separates the two groups?" — a strictly larger
# question, answered from the data, and therefore Family B forever.
#
# Correlated with nothing in Family A, and never pooled with it: family_adjust()
# keeps the two families apart precisely so a discovery result cannot spend the
# confirmatory panel's power.

#' Per-sample "altered at this peak", as a logical matrix.
#'
#' ⚠️ all_lesions.conf_*.txt is UNSIGNED: a deletion peak's thresholded call is 1 or 2,
#' not -1 or -2 — the peak's direction lives in its Unique Name, i.e. peaks$direction,
#' never in the call's sign. An earlier version took direction from the sign of a
#' group-mean difference and tested `sign * call >= 1`, which on this matrix is never
#' true for a negative sign: every locus MORE frequent in the second group scored FALSE
#' for everyone, so the data-driven panel could only ever find loci enriched in the
#' first group, and its "amp"/"del" column was a group label under the wrong name. The
#' synthetic test matrices were signed, so nothing failed.
#'
#' `signed` is DECLARED, never inferred from the values: a signed matrix in which no
#' sample happens to carry a negative call is indistinguishable from an unsigned one.
#' Unsigned (the default, all_lesions): altered is call >= 1, whatever the direction.
#' Signed with `peak_dir`: only concordant calls count (+ at amp, - at del).
.peak_altered <- function(mat, peak_dir = NULL, signed = FALSE) {
  if (!signed || is.null(peak_dir)) return(abs(mat) >= 1)
  sweep(mat, 2, ifelse(peak_dir == "del", -1, 1), `*`) >= 1
}

#' Choose the n_select most group-discriminating peaks from an ALTERED matrix.
#'
#' DELIBERATELY label-dependent — the reason the result can never be treated as
#' pre-specified. Criterion is the difference in alteration FREQUENCY between the two
#' groups: magnitude ranks the peaks, sign records which group each one is enriched in
#' (`enrich` +1 = first level, -1 = second). A peak's amp/del direction is a property of
#' the peak, not of this selection.
.select_on_altered <- function(A, grp, n_select) {
  lv <- levels(droplevels(factor(grp)))
  if (length(lv) != 2L) stop("select_panel_loci(): need exactly 2 groups, got ", length(lv))
  d <- colMeans(A[grp == lv[1], , drop = FALSE], na.rm = TRUE) -
       colMeans(A[grp == lv[2], , drop = FALSE], na.rm = TRUE)
  d[!is.finite(d)] <- 0
  # A zero-difference peak is enriched in neither group, so it is not selectable.
  cand <- which(d != 0)
  if (!length(cand)) return(list(idx = integer(0), enrich = numeric(0)))
  idx <- cand[order(abs(d[cand]), decreasing = TRUE)][seq_len(min(n_select, length(cand)))]
  list(idx = idx, enrich = sign(d[idx]))
}

select_panel_loci <- function(mat, grp, n_select = attend_scna$select_n, peak_dir = NULL) {
  .select_on_altered(.peak_altered(mat, peak_dir), grp, n_select)
}

#' Per-sample data-driven panel score: the fraction of selected loci at which the
#' sample looks like the FIRST group — altered where that group is enriched, unaltered
#' where the second group is. Both kinds of locus contribute, so the group-mean
#' difference equals the mean |frequency difference| over the selected loci.
.selected_panel_score <- function(A, sel) {
  if (!length(sel$idx)) return(rep(NA_real_, nrow(A)))
  sub <- A[, sel$idx, drop = FALSE]
  rowMeans(sweep(sub, 2, sel$enrich > 0, `==`), na.rm = TRUE)
}

#' Group-mean difference in data-driven panel score, for one labelling.
.selected_panel_stat <- function(A, grp, n_select) {
  lv  <- levels(droplevels(factor(grp)))
  sel <- .select_on_altered(A, grp, n_select)
  if (!length(sel$idx)) return(0)
  s   <- .selected_panel_score(A, sel)
  mean(s[grp == lv[1]], na.rm = TRUE) - mean(s[grp == lv[2]], na.rm = TRUE)
}

#' Permutation test for a data-driven panel, with the selection re-run inside every
#' replicate.
#'
#' ⚠️ THE VALIDITY IS ENTIRELY IN WHERE THE SELECTION HAPPENS. The selection is called
#' INSIDE the replicate, on the PERMUTED labels. Hoisting it out — select once, then
#' permute the scores — is the invalid "double dipping" version (Kriegeskorte et al.,
#' Nat Neurosci 2009): the panel would be chosen with the real labels, so the permuted
#' scores could never reproduce the optimism the selection introduced, the null would be
#' far too narrow, and p would be near zero on pure noise. That hoist looks like an
#' obvious speed-up and is invisible in a diff. test_nested_selection_permutation.R pins
#' the CALIBRATION (uniform p under the null), which is what actually catches it.
#'
#' The altered matrix is computed ONCE, outside the replicates: it uses no labels.
#'
#' ONE-SIDED by construction: the statistic is a mean |frequency difference| and
#' cannot be negative, so comparing |null| to |obs| would test a hypothesis it cannot
#' express.
perm_test_selected_panel <- function(mat, grp,
                                     n_select = attend_scna$select_n,
                                     B        = attend_scna$select_perm_B,
                                     seed     = 1,
                                     peak_dir = NULL) {
  keep <- !is.na(grp)
  A    <- .peak_altered(mat, peak_dir)[keep, , drop = FALSE]
  grp  <- droplevels(factor(grp[keep]))
  if (nlevels(grp) != 2L)
    stop("perm_test_selected_panel(): need exactly 2 groups, got ", nlevels(grp))

  obs <- .selected_panel_stat(A, grp, n_select)
  set.seed(seed)
  null <- replicate(B, .selected_panel_stat(A, sample(grp), n_select))

  sel <- .select_on_altered(A, grp, n_select)
  # Leave-one-out selection stability. A panel whose membership turns over when one
  # patient is dropped is not a finding, and at n ~ 9 per group that is the likely
  # outcome — so it is reported beside p rather than left for a reader to wonder about.
  loo <- table(unlist(lapply(seq_len(nrow(A)), function(i)
    colnames(A)[.select_on_altered(A[-i, , drop = FALSE], grp[-i], n_select)$idx])))

  loci <- colnames(A)[sel$idx]
  list(stat        = obs,
       p           = (1 + sum(null >= obs)) / (B + 1),
       B           = B,
       loci        = loci,
       peak_dir    = if (is.null(peak_dir)) rep(NA_character_, length(loci))
                     else unname(peak_dir[sel$idx]),
       enriched_in = levels(grp)[ifelse(sel$enrich > 0, 1L, 2L)],
       score       = stats::setNames(.selected_panel_score(A, sel), rownames(A)),
       loo_frac    = as.numeric(loo[loci]) / nrow(A))
}

#' Global two-group test over the WHOLE peak set, selecting nothing.
#'
#' The complement to both panels: total burden per sample (any call of |value| >= 1 at
#' any peak), through the fixed-panel permutation. No selection, so no selection to
#' correct for, and it stays answerable when both panel routes come back null — it asks
#' "is the profile different at all?", which neither a serous-aimed panel nor a
#' discriminating-peak panel can answer.
perm_test_global_burden <- function(mat, grp, B = attend_scna$perm_B, seed = 1) {
  burden <- rowMeans(abs(mat) >= 1, na.rm = TRUE)
  perm_test_two_group(burden, grp, B = B, seed = seed)
}

#' Burden control for the pre-specified panel: is the group difference at the 12
#' serous loci larger than at RANDOM loci of the same make-up?
#'
#' perm_test_two_group() cannot separate "serous-like" from "more aneuploid". Here the
#' observed MMRp-high minus MMRd-high difference in panel score is set against the same
#' difference over B random panels drawn from the pooled run's peaks with the SAME
#' number of amplification and deletion loci the panel resolved. A burden difference
#' moves every random panel too, so it moves the null, not the contrast.
#'
#' The resampling unit is the LOCUS, not the patient, so p here says "specific to these
#' loci", not "a patient-level effect" — it is the sensitivity analysis beside the
#' primary endpoint, never a replacement for it.
panel_specificity_test <- function(mat, peaks, grp, panel = attend_scna$panel,
                                   B = attend_scna$perm_B, seed = 1) {
  mp <- match_panel_peaks(peaks, panel)
  mp <- mp[!is.na(mp$peak_id) & mp$peak_id %in% colnames(mat), , drop = FALSE]
  keep <- !is.na(grp)
  g <- droplevels(factor(grp[keep]))
  if (!nrow(mp) || nlevels(g) != 2L) return(NULL)

  pdir <- peaks$direction[match(colnames(mat), peaks$peak_id)]
  A    <- .peak_altered(mat, pdir)[keep, , drop = FALSE]
  lv   <- levels(g)
  diff_of <- function(ids) {
    s <- rowMeans(A[, ids, drop = FALSE], na.rm = TRUE)
    mean(s[g == lv[1]]) - mean(s[g == lv[2]])
  }
  n_amp <- sum(mp$direction == "amp"); n_del <- sum(mp$direction == "del")
  pool_amp <- colnames(A)[pdir %in% "amp"]; pool_del <- colnames(A)[pdir %in% "del"]
  if (length(pool_amp) < n_amp || length(pool_del) < n_del) return(NULL)

  obs <- diff_of(mp$peak_id)
  set.seed(seed)
  null <- replicate(B, diff_of(c(sample(pool_amp, n_amp), sample(pool_del, n_del))))
  ctr  <- mean(null)
  list(stat = obs, null_mean = ctr,
       null_lo = unname(stats::quantile(null, 0.025)),
       null_hi = unname(stats::quantile(null, 0.975)),
       p = (1 + sum(abs(null - ctr) >= abs(obs - ctr))) / (B + 1),
       B = B, n_amp = n_amp, n_del = n_del)
}

# --- the two questions, as tables ----------------------------------------------

#' Question 1 — which loci recur in one stratum, and how firmly.
#'
#' Candidate loci are the UNION of the pooled run's peaks and the stratum run's own
#' peaks, matched by wide-limit overlap on chromosome and direction. The pooled run
#' alone would bury a peak private to a nine-patient stratum under its genome-wide
#' threshold, which is the reason the stratum run exists.
#'
#' "Recurrent" is GISTIC's own call on the STRATUM run: q <= fdr there. A frequency
#' cut ("altered in half the tumours") is not recurrence — aneuploidy-high tumours carry
#' so many arm-level events that a high frequency at a pooled peak is what burden alone
#' predicts. Frequency, its Wilson CI and leave-one-out retention are reported beside
#' the call, as its support, never in place of it.
#'
#' Per-sample calls come from the POOLED run whenever the locus has a pooled peak
#' (one background model for every comparison that follows). Only a locus with no
#' pooled peak falls back to the stratum run's own calls, and `calls_from` says so —
#' those rows are within-stratum only and cannot enter a between-group comparison.
q1_recurrence_table <- function(pooled, stratum, grp, target = "MMRd-high",
                                fdr = 0.25, loo = NULL, overlap = peaks_overlap) {
  if (is.null(pooled) || is.null(stratum)) return(NULL)
  pk <- pooled$peaks; sk <- stratum$peaks

  # Stratum peak -> the most significant overlapping pooled peak, if any.
  map <- vapply(seq_len(nrow(sk)), function(i) {
    ok <- vapply(seq_len(nrow(pk)), function(j) isTRUE(overlap(sk[i, ], pk[j, ])), logical(1))
    if (!any(ok)) NA_character_ else pk$peak_id[ok][which.min(pk$q_value[ok])]
  }, character(1))

  ft  <- scna_freq_table(pooled$mat, grp)
  ft  <- ft[as.character(ft$scna_group) == target, , drop = FALSE]
  fs  <- scna_freq_table(stratum$mat, factor(rep(target, nrow(stratum$mat))))
  row_for <- function(tbl, id) tbl[match(id, tbl$peak_id), c("n", "n_altered", "freq", "ci_lo", "ci_hi")]

  from_stratum <- data.frame(
    locus = sk$descriptor, direction = sk$direction,
    found_in = ifelse(is.na(map), paste(target, "run only"), "both runs"),
    q_stratum = sk$q_value, recurrent = !is.na(sk$q_value) & sk$q_value <= fdr,
    stratum_peak_id = sk$peak_id, pooled_peak_id = map, stringsAsFactors = FALSE)
  from_stratum <- cbind(from_stratum,
                        ifelse_rows(is.na(map), row_for(fs, sk$peak_id), row_for(ft, map)))
  from_stratum$calls_from <- ifelse(is.na(map), paste(target, "run"), "pooled run")

  rest <- pk[!pk$peak_id %in% map, , drop = FALSE]
  from_pooled <- data.frame(
    locus = rest$descriptor, direction = rest$direction, found_in = "pooled run only",
    q_stratum = NA_real_, recurrent = FALSE, stratum_peak_id = NA_character_,
    pooled_peak_id = rest$peak_id, stringsAsFactors = FALSE)
  from_pooled <- cbind(from_pooled, row_for(ft, rest$peak_id))
  from_pooled$calls_from <- if (nrow(rest)) "pooled run" else character(0)

  out <- rbind(from_stratum, from_pooled)
  out$loo_retained_frac <- if (is.null(loo)) NA_real_
                           else loo$retained_frac[match(out$stratum_peak_id, loo$peak_id)]
  out <- out[order(!out$recurrent, out$q_stratum, -out$freq, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Row-wise choice between two equally-shaped data frames.
ifelse_rows <- function(test, yes, no) {
  out <- no; out[test, ] <- yes[test, ]; rownames(out) <- NULL; out
}

#' Question 2 — the loci Q1 called recurrent, scored in the reference stratum.
#'
#' Restricted to Q1's recurrent loci because the question is conditional ("if so, are
#' THEY the same"). Both frequencies come from the pooled run's calls. A locus with no
#' pooled peak is kept with NA and a reason: it cannot be scored in the reference group
#' on the same background model, and dropping it would hide the stratum's most specific
#' finding. No verdict column — at n ~ 9 the intervals are the answer.
q2_same_loci_table <- function(q1, pooled, grp, target = "MMRd-high", ref = "MMRp-high") {
  if (is.null(q1)) return(NULL)
  r <- q1[q1$recurrent, , drop = FALSE]
  if (!nrow(r)) return(r[, 0, drop = FALSE])
  ft <- scna_freq_table(pooled$mat, grp)
  ft <- ft[as.character(ft$scna_group) == ref, , drop = FALSE]
  m  <- match(r$pooled_peak_id, ft$peak_id)
  data.frame(locus = r$locus, direction = r$direction, q_stratum = r$q_stratum,
             freq_target = r$freq, ci_lo_target = r$ci_lo, ci_hi_target = r$ci_hi,
             freq_ref = ft$freq[m], ci_lo_ref = ft$ci_lo[m], ci_hi_ref = ft$ci_hi[m],
             n_ref = ft$n[m], diff = r$freq - ft$freq[m],
             note = ifelse(is.na(r$pooled_peak_id), "no pooled-run peak: not scorable in ref",
                           ifelse(r$calls_from == "pooled run", "", r$calls_from)),
             stringsAsFactors = FALSE)
}

#' How concordant would two groups be if they were the SAME population, at these n?
#'
#' Spearman's rho between per-peak frequencies has no meaningful fixed benchmark here:
#' peaks near 0% in both groups, shared arm-level events and correlated peaks on one arm
#' all push it up, while n ~ 9 pulls it down. So it is read against two references:
#'  - `null_*`: rho between two random splits of the SAME patients (target + ref pooled,
#'    labels shuffled, group sizes kept). This is the ceiling "the same" reaches at this
#'    n. Permuting PATIENTS keeps every peak-peak correlation, so arm structure is in the
#'    null too. An observed rho inside this range is consistent with "the same"; below
#'    its 2.5% quantile it is less concordant than chance splits of one population.
#'  - `ref_rho`: rho of target against a group expected to differ (default MMRp-low), the
#'    floor.
#' `cols` restricts to a peak subset (e.g. Q1's recurrent loci). ⚠️ Those loci were
#' chosen for being frequent in the target, so the target's frequencies there are biased
#' upward and the reference's regress; that bias pushes rho DOWN, toward "different", and
#' the shuffled splits do not carry it.
concordance_null <- function(mat, grp, target = "MMRd-high", ref = "MMRp-high",
                             floor_group = "MMRp-low", cols = NULL, B = 1000, seed = 1) {
  A <- abs(mat) >= 1
  if (!is.null(cols)) A <- A[, intersect(cols, colnames(A)), drop = FALSE]
  i1 <- which(grp == target); i2 <- which(grp == ref)
  if (ncol(A) < 3 || !length(i1) || !length(i2)) return(NULL)
  f   <- function(rows) colMeans(A[rows, , drop = FALSE], na.rm = TRUE)
  rho <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))

  obs  <- rho(f(i1), f(i2))
  pool <- c(i1, i2); k <- length(i1)
  set.seed(seed)
  null <- replicate(B, { s <- sample(pool); rho(f(s[seq_len(k)]), f(s[-seq_len(k)])) })
  i3 <- which(grp == floor_group)
  list(rho = obs, n_peaks = ncol(A),
       null_median = stats::median(null, na.rm = TRUE),
       null_lo = unname(stats::quantile(null, 0.025, na.rm = TRUE)),
       null_hi = unname(stats::quantile(null, 0.975, na.rm = TRUE)),
       frac_null_below = mean(null <= obs, na.rm = TRUE),
       floor_group = floor_group,
       ref_rho = if (length(i3)) rho(f(i1), f(i3)) else NA_real_,
       B = B)
}
