# relabel_clusters_by_burden() — the paper's x-axis order.
#
# WHY THIS EXISTS. Kandoth et al. Fig. 1a runs left to right from copy-number quiet to
# copy-number high. Measured on TCGA's own published CNA_CLUSTER_K4 over their hg19 segments,
# mean |log2| per cluster is 0.0043 / 0.052 / 0.099 / 0.216 — monotonic over a 50x range, and
# their "Copy-number high (Serous-like)" subtype is 60/60 inside cluster 4. TCGA numbered by
# burden. cutree() numbers by order of first appearance in the dendrogram, which is an
# artefact of merge order, so a faithful clustering still drew its blocks in a meaningless
# sequence. On the real TCGA peak matrix this relabelling lifts diagonal agreement with the
# published labels from 23.2% to 41.1% without touching the partition (ARI is unchanged at
# 0.222, being label-invariant) — it is presentation, and it is the presentation the paper has.
#
# Base R only, so this runs in the bootstrap env.
source(file.path("code", "attend_plots.R"))

stopifnot(is.function(relabel_clusters_by_burden))

## --- ascending order, and the map that goes with it -------------------------
## Three samples per cluster. Cluster "a" is the loudest, "c" the quietest, so the expected
## relabelling is c -> 1, b -> 2, a -> 3.
mk <- function(v, n) matrix(v, nrow = n, ncol = 4, byrow = TRUE)
m <- rbind(mk(c( 2, -2,  2, -2), 3),   # a: mean |v| = 2
           mk(c( 1,  0, -1,  0), 3),   # b: mean |v| = 0.5
           mk(c( 0,  0,  0,  0), 3))   # c: mean |v| = 0  (flat)
rownames(m) <- paste0("s", 1:9)
cl <- stats::setNames(rep(c("a", "b", "c"), each = 3), rownames(m))

r <- relabel_clusters_by_burden(cl, m)
stopifnot(identical(unname(r$map[c("c", "b", "a")]), c("1", "2", "3")))
stopifnot(identical(as.character(r$cluster[c("s1", "s4", "s7")]), c("3", "2", "1")))
## burden is reported against the OLD labels, so a caller can print what drove the order
stopifnot(abs(r$burden[["a"]] - 2) < 1e-9,
          abs(r$burden[["b"]] - 0.5) < 1e-9,
          abs(r$burden[["c"]] - 0)   < 1e-9)
## and the relabelled clusters really do ascend
b <- tapply(rowMeans(abs(m)), r$cluster, mean)
stopifnot(!is.unsorted(b))

## levels are the full 1..k sequence, so a heatmap column split has a stable order even when
## a cluster is empty in some subset of the data.
stopifnot(is.factor(r$cluster), identical(levels(r$cluster), c("1", "2", "3")))

## --- the partition is PRESERVED: only labels move ---------------------------
## Two samples share a new label if and only if they shared an old one. If this ever fails the
## relabelling is not a relabelling, it is a re-clustering.
# unname(): outer() on a NAMED vector carries dimnames while outer() on a bare character
# vector does not, so identical() would compare attributes rather than the partition.
same_old <- unname(outer(cl, cl, `==`))
same_new <- unname(outer(as.character(r$cluster), as.character(r$cluster), `==`))
stopifnot(identical(dim(same_old), dim(same_new)), all(same_old == same_new))

## --- a cluster with no row in the matrix is ordered LAST, never slot 1 ------
## Otherwise an unmeasured cluster would silently claim "quietest", which is a claim about
## data that does not exist.
cl2 <- c(cl, s99 = "ghost")
r2  <- relabel_clusters_by_burden(cl2, m)          # s99 is absent from m
stopifnot(identical(unname(r2$map[["ghost"]]), "4"),
          is.na(r2$burden[["ghost"]]))

## --- degenerate inputs are no-ops, not errors -------------------------------
stopifnot(identical(relabel_clusters_by_burden(cl, NULL)$cluster, cl))
e <- relabel_clusters_by_burden(stats::setNames(character(0), character(0)), m)
stopifnot(length(e$cluster) == 0L)
## an unnamed cluster vector cannot be matched to matrix rows — return it untouched rather
## than guessing a row order, which would assign burden to the wrong samples.
stopifnot(identical(relabel_clusters_by_burden(unname(cl), m)$cluster, unname(cl)))

## --- ties keep a stable, deterministic order --------------------------------
## Two clusters with identical burden must not reorder between runs, or the same data would
## produce two different figures.
m3 <- rbind(mk(c(1, -1, 1, -1), 2), mk(c(1, -1, 1, -1), 2))
rownames(m3) <- paste0("t", 1:4)
cl3 <- stats::setNames(rep(c("x", "y"), each = 2), rownames(m3))
stopifnot(identical(relabel_clusters_by_burden(cl3, m3)$map,
                    relabel_clusters_by_burden(cl3, m3)$map))

cat("relabel_clusters_by_burden: ALL PASS\n")
