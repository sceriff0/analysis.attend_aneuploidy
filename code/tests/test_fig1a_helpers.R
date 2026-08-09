# Base-R unit tests (no tidyverse) — runnable in the bootstrap env.
source(file.path("code", "attend_plots.R"))  # base-R-sourceable: no top-level library()

## --- feature_positions: genome bins ---
fp <- feature_positions(c("chr2:1000000", "chr1:0", "chr1:1000000"))
stopifnot(is.data.frame(fp), nrow(fp) == 3)
# ordered by (chrom, coord): chr1:0, chr1:1e6, chr2:1e6
stopifnot(identical(as.character(fp$feature), c("chr1:0","chr1:1000000","chr2:1000000")))
stopifnot(identical(as.character(fp$chrom), c("chr1","chr1","chr2")))
stopifnot(fp$coord == c(0, 1000000, 1000000))

## --- feature_positions: arms ---
fa <- feature_positions(c("17q", "1p", "1q"))
stopifnot(identical(as.character(fa$feature), c("1p","1q","17q")))     # chr1 before chr17, p before q
stopifnot(identical(as.character(fa$chrom), c("chr1","chr1","chr17")))

## --- feature_positions: GISTIC genes via cytoband ---
cb <- c(MYC = "8q24.21", TP53 = "17p13.1", CCNE1 = "19q12")
fg <- feature_positions(c("CCNE1","MYC","TP53"), cytoband = cb)
stopifnot(identical(as.character(fg$chrom), c("chr8","chr17","chr19")))  # genome order 8,17,19
stopifnot(all(c("MYC","TP53","CCNE1") %in% fg$feature))

cat("feature_positions: ALL PASS\n")

## --- fig1a_heatmap: smoke (skips cleanly without ComplexHeatmap) ---
set.seed(1)
mat <- matrix(sample(c(-2,-1,0,1,2), 6*4, replace = TRUE), nrow = 4,
              dimnames = list(paste0("S", 1:4), c("chr1:0","chr1:1000000","chr8:1000000","17q","19q12","2p")))
fp  <- feature_positions(colnames(mat), cytoband = c("19q12"="19q12"))
cls <- c(S1=1, S2=1, S3=2, S4=2)
res <- fig1a_heatmap(mat, fp, cls, ann_col = NULL, value_type = "thresholded",
                     main = "smoke")   # must not error whether or not CH is installed
stopifnot(is.null(res) || inherits(res, "Heatmap") || inherits(res, "HeatmapList"))
cat("fig1a_heatmap: SMOKE OK\n")

## --- km_by_cluster: computes a log-rank p from cluster + PFS (survival is present) ---
# attend_plots.R does not define recode_event (it lives in attend_classes.R, which
# can't be sourced here); provide a base-R stub matching its text->0/1 contract.
if (!exists("recode_event"))
  recode_event <- function(x) as.integer(grepl("pd|death|dead|progress|deceased",
                                               tolower(as.character(x))))
assignment <- data.frame(ID = paste0("bc", 1:8), cnv_cluster = rep(1:2, each = 4), stringsAsFactors = FALSE)
cw     <- data.frame(barcode = paste0("bc", 1:8), pid = paste0("p", 1:8), stringsAsFactors = FALSE)
master <- data.frame(pid = paste0("p", 1:8),
                     gianlu__PFS_MONTHS = c(5,6,7,8, 20,22,25,30),
                     gianlu__PFS_EVENT  = c("PD","PD","PD","Death", "No event","No event","PD","No event"),
                     stringsAsFactors = FALSE)
km <- km_by_cluster(assignment, cw, master,
                    time_col = "gianlu__PFS_MONTHS", event_col = "gianlu__PFS_EVENT")
stopifnot(is.null(km) || (is.list(km) && !is.null(km$p_logrank) && is.numeric(km$p_logrank)))
cat("km_by_cluster: OK\n")
