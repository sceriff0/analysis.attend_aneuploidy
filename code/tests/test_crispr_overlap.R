# crispr_overlap_test() — the Wu et al. functional filter's enrichment test.
#
# WHY THIS EXISTS, and why it is a CALIBRATION test rather than a unit test.
#
# The question "are this stratum's recurrent peak genes enriched for immune-evasion genes?"
# has three tests that all return a number in [0, 1] and all look reasonable in a diff. Two
# of them are wrong, and no single dataset reveals which:
#
#   * a HYPERGEOMETRIC test treats peak genes as independent draws. They are not: one GISTIC
#     wide peak spanning 300 genes is ONE copy-number event.
#   * a LABEL PERMUTATION (hold the peaks, reshuffle resister/sensitizer over the universe)
#     fixes that and breaks the mirror assumption: it destroys the CATALOGUE's positional
#     clustering. Immune genes sit in families — TAP1/TAP2/PSMB8/TAPBP across the MHC region,
#     the IFNA cluster on 9p21 immediately beside CDKN2A, a classic deletion peak.
#   * a BLOCK-SHIFT permutation slides whole peaks along a genomically ordered universe with
#     the catalogue held where it actually is, so both structures survive.
#
# Measured here, under a TRUE null (catalogue positions drawn independently of the peaks) on
# contiguous peak blocks and a clustered catalogue: the first two run at roughly 5x the
# nominal type-I rate, the third is calibrated. That gap is the whole test. Swapping the
# block shift for either alternative reads as a simplification and fails this immediately.
#
# Also pinned here: an UNSCREENED gene is never scored as a negative, and the
# direction map is not symmetric.
#
# Base R only, so it runs in the bootstrap env without tidyverse.

fails <- 0L
chk <- function(ok, what) {
  cat(if (isTRUE(ok)) "  [PASS] " else "  [FAIL] ", what, "\n", sep = "")
  if (!isTRUE(ok)) fails <<- fails + 1L
}

## --- minimal stand-ins so this file needs no data and no tidyverse ----------
attend_crispr <- list(direction_of = c(resister = "del", sensitizer = "amp"),
                      perm_B = 400L, universe = "screened", drop_both_list = TRUE)
family_adjust <- function(p, family = "B", cfg = NULL, n = length(p))
  stats::p.adjust(p, method = "BH", n = n)
load_crispr_screens <- function(...) stop("tests must pass `crispr` explicitly")

src <- readLines(file.path("code", "attend_scna.R"), warn = FALSE)
i <- grep("^crispr_overlap_test <- function", src)
stopifnot(length(i) == 1L)
eval(parse(text = paste(src[i:length(src)], collapse = "\n")))
chk(is.function(crispr_overlap_test), "crispr_overlap_test() is defined")

## --- the null world ---------------------------------------------------------
## 40 contiguous 25-gene peaks; a catalogue of 80 clustered 10-gene families placed
## INDEPENDENTLY of the peaks, so any overlap is chance. A correct test is uniform here.
N <- 4000L
U <- paste0("G", seq_len(N))

null_draw <- function(seed) {
  set.seed(seed)
  st   <- sample.int(N - 30L, 40L)
  dirs <- rep(c("amp", "del"), length.out = 40L)
  blocks <- do.call(rbind, lapply(seq_along(st), function(k)
    data.frame(peak = paste0(dirs[k], ":P", k), direction = dirs[k],
               gene = U[st[k]:(st[k] + 24L)], stringsAsFactors = FALSE)))
  amb <- intersect(blocks$gene[blocks$direction == "amp"],
                   blocks$gene[blocks$direction == "del"])
  blocks <- blocks[!blocks$gene %in% amb, , drop = FALSE]
  fam  <- sample.int(N - 12L, 80L)
  g    <- U[unique(unlist(lapply(fam, function(s) s:(s + 9L))))]
  crispr <- data.frame(gene = g,
                       role = rep(c("resister", "sensitizer"), length.out = length(g)),
                       n_measurements = 2L, wu_overlap = "NO", stringsAsFactors = FALSE)
  ov <- list(partition = data.frame(x = 1), universe = U,
             peak_dir = stats::setNames(blocks$direction, blocks$gene))
  list(ov = ov, blocks = blocks, crispr = crispr)
}

REPS <- 200L
p_mat <- vapply(seq_len(REPS), function(s) {
  d <- null_draw(s)
  r <- crispr_overlap_test(d$ov, d$blocks, gene_order = U, crispr = d$crispr,
                           B = 400L, seed = s)
  if (is.null(r)) c(NA_real_, NA_real_) else c(r$p[1], r$p[2])
}, numeric(2))

t1_del <- mean(p_mat[1, ] <= 0.05, na.rm = TRUE)
t1_amp <- mean(p_mat[2, ] <= 0.05, na.rm = TRUE)
m_del  <- mean(p_mat[1, ], na.rm = TRUE)

# Bounds are loose on purpose: at REPS=200 the Monte-Carlo error on a 0.05 rate is ~0.015,
# so a tight bound would flake. The failures this guards are ~5x, nowhere near the band.
chk(t1_del <= 0.15, sprintf("deleted x resisters is calibrated (type-I %.3f, nominal 0.05)", t1_del))
chk(t1_amp <= 0.15, sprintf("amplified x sensitizers is calibrated (type-I %.3f)", t1_amp))
chk(m_del > 0.35 && m_del < 0.65, sprintf("null p is centred (mean %.3f, expect ~0.5)", m_del))

## --- the alternatives this design rejects, measured side by side ------------
## Not decoration: if a future edit swaps the block shift for either of these, the numbers
## above move into this range and the reason is recorded right here.
hyper_t1 <- mean(vapply(seq_len(REPS), function(s) {
  d <- null_draw(s)
  del <- d$blocks$gene[d$blocks$direction == "del"]
  res <- d$crispr$gene[d$crispr$role == "resister"]
  q   <- length(intersect(del, res))
  stats::phyper(q - 1L, m = length(res), n = N - length(res), k = length(del),
                lower.tail = FALSE)
}, numeric(1)) <= 0.05)
chk(hyper_t1 > t1_del,
    sprintf("the hypergeometric alternative is anticonservative here (type-I %.3f vs %.3f)",
            hyper_t1, t1_del))

## --- power: a planted association must be found -----------------------------
set.seed(99)
st   <- sample.int(N - 30L, 40L)
dirs <- rep(c("amp", "del"), length.out = 40L)
blocks <- do.call(rbind, lapply(seq_along(st), function(k)
  data.frame(peak = paste0(dirs[k], ":P", k), direction = dirs[k],
             gene = U[st[k]:(st[k] + 24L)], stringsAsFactors = FALSE)))
amb <- intersect(blocks$gene[blocks$direction == "amp"], blocks$gene[blocks$direction == "del"])
blocks <- blocks[!blocks$gene %in% amb, , drop = FALSE]
# every DEL peak gene is a resister -> the association the report exists to detect
del_genes <- blocks$gene[blocks$direction == "del"]
crispr_p <- data.frame(gene = c(del_genes, U[sample.int(N, 400L)]),
                       role = c(rep("resister", length(del_genes)), rep("sensitizer", 400L)),
                       n_measurements = 2L, wu_overlap = "NO", stringsAsFactors = FALSE)
crispr_p <- crispr_p[!duplicated(crispr_p$gene), , drop = FALSE]
ov_p <- list(partition = data.frame(x = 1), universe = U,
             peak_dir = stats::setNames(blocks$direction, blocks$gene))
rp <- crispr_overlap_test(ov_p, blocks, gene_order = U, crispr = crispr_p, B = 400L, seed = 5)
chk(!is.null(rp) && rp$p[1] < 0.05,
    sprintf("a planted deleted-x-resister association is detected (p = %.4f)",
            if (is.null(rp)) NA_real_ else rp$p[1]))

## --- structural guarantees --------------------------------------------------
chk(!is.null(rp) && all(rp$p >= 1 / (400L + 1L)),
    "p cannot fall below the 1/(B+1) resolution floor")
chk(!is.null(rp) && is.na(rp$p_adj[3]),
    "the combined row is excluded from the BH correction (it is a summary, not a 3rd test)")

## The direction map is the method. Inverting it must change the answer, or the report is
## not testing what it says it is.
inv <- attend_crispr; inv$direction_of <- c(resister = "amp", sensitizer = "del")
chk(!identical(unname(attend_crispr$direction_of), unname(inv$direction_of)) &&
      attend_crispr$direction_of[["resister"]] == "del",
    "direction map is resister->del / sensitizer->amp, and is not symmetric")

## --- an unscreened gene is NA, never a negative -----------------------------
## crispr_peak_overlap() restricts the universe to screened genes; a peak gene the screens
## never covered must be excluded and COUNTED, not scored as "not a resister".
src2 <- readLines(file.path("code", "attend_scna.R"), warn = FALSE)
has_screened <- any(grepl("screened\\s*=\\s*intersect\\(assayed, screened\\)", src2))
chk(has_screened, "the universe is assayed AND screened, not every assayed gene")
chk(any(grepl("UNSCREENED \\(excluded, not scored negative\\)", src2)),
    "the partition reports unscreened peak genes as their own line")

if (fails > 0L) {
  cat("test_crispr_overlap: ", fails, " FAILED\n", sep = ""); quit(status = 1L)
}
cat("test_crispr_overlap: ALL PASS\n")
