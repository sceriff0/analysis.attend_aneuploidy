# load_gistic_continuous_at() — the Fig-1a DISPLAY matrix, which is NOT the clustering matrix.
#
# WHY THIS TEST EXISTS. Kandoth et al. use two different matrices and report 09 conflated them:
# it clustered the peak-restricted thresholded calls (correct, Suppl. Methods S2) and then DREW
# that same matrix, so the heatmap was a few hundred peak genes on a five-step colour scale
# instead of the paper's genome-wide continuous landscape. The two properties that make this
# loader the display matrix are exactly the two an "optimisation" would undo — no peak
# restriction, and no integer coercion — so both are pinned here.
#
# NOTE: load_wes_results.R does `library(tidyverse)` at top, so this test requires tidyverse.
# In the bootstrap/dev env where tidyverse is absent, verify instead with
#   Rscript -e 'invisible(parse("code/load_wes_results.R")); cat("PARSE OK\n")'
# and let the cluster run the real check.
source(file.path("code", "load_wes_results.R"))

dir.create(td <- file.path(tempdir(), "gistic_cont"), showWarnings = FALSE)

# all_data_by_genes.txt: col1 Gene Symbol, col2 Gene ID, col3 Cytoband, then samples.
# CONTINUOUS values, deliberately non-integer, and spanning genes that are NOT in any peak.
dat <- data.frame(
  `Gene Symbol` = c("MYC", "TP53", "PTEN", "SOMEGENE"),
  `Gene ID`     = c(1, 2, 3, 4),
  Cytoband      = c("8q24.21", "17p13.1", "10q23.31", "4q13.2"),
  S1 = c( 0.8412, -0.9105, -0.4400,  0.0123),
  S2 = c( 0.1200,  0.0050, -1.2345, -0.0500),
  check.names = FALSE)
write.table(dat, file.path(td, "all_data_by_genes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# A peak-gene list that covers ONLY MYC. The thresholded loader restricts to it; this one
# must ignore it entirely, or the display matrix loses the genome and becomes the clustering
# matrix again — the exact regression this loader was written to end.
amp <- data.frame(`amp:8q24.21` = c("cytoband", "q value", "MYC"), check.names = FALSE)
write.table(amp, file.path(td, "amp_genes.conf_99.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cnv <- list(gistic = list(dir = basename(td),
                          all_data_glob        = "*all_data_by_genes*.txt",
                          all_thresholded_glob = "*all_thresholded*.txt",
                          amp_genes_glob       = "*amp_genes*.txt",
                          del_genes_glob       = "*del_genes*.txt"),
            seg = list(id_strip = "-1TAD104|_tumor_only"))

x <- tryCatch(load_gistic_continuous_at(td, cnv),
              error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
stopifnot(!is.null(x))

## samples are ROWS (ID column), genes are columns — the shape cluster_arm_matrix() and
## fig1a_heatmap() both expect.
stopifnot("ID" %in% names(x), nrow(x) == 2L, all(c("S1", "S2") %in% x$ID))

## NOT peak-restricted: all four genes survive even though only MYC is in a peak.
genes <- setdiff(names(x), "ID")
stopifnot(length(genes) == 4L,
          all(c("MYC", "TP53", "PTEN", "SOMEGENE") %in% genes))

## CONTINUOUS: values are doubles and keep their decimals. as.integer() would have turned
## 0.8412 into 0 and -0.9105 into 0, silently rebuilding a near-empty discrete matrix.
v <- x[[which(x$ID == "S1"), "MYC"]]
stopifnot(is.numeric(v), !is.integer(v), abs(v - 0.8412) < 1e-9)
stopifnot(abs(x[[which(x$ID == "S2"), "PTEN"]] - (-1.2345)) < 1e-9)
## and nothing rounded to zero
stopifnot(all(vapply(x[genes], function(col) any(col %% 1 != 0), logical(1))))

## feature_pos carries the cytoband per gene, keyed by the SAME make.unique() names as the
## columns, so feature_positions(colnames(mat), cytoband = ...) resolves every row.
fp <- attr(x, "feature_pos")
stopifnot(!is.null(fp), all(c("feature", "cytoband") %in% names(fp)),
          setequal(fp$feature, genes),
          fp$cytoband[fp$feature == "PTEN"] == "10q23.31")

## End-to-end with the figure layer: every displayed gene must land on the genome, in order.
## This is the join that silently drops every row if the cytoband keys come from the wrong
## GISTIC file — the reason report 09 reads feature_pos off the DISPLAY matrix.
source(file.path("code", "attend_plots.R"))
pos <- feature_positions(genes, cytoband = setNames(fp$cytoband, fp$feature))
stopifnot(nrow(pos) == 4L,
          all(as.character(pos$chrom) %in% paste0("chr", c(4, 8, 10, 17))),
          !any(is.na(pos$coord)),
          identical(as.character(pos$chrom), as.character(sort(pos$chrom))))

## A missing file is knit-safe: NULL, so the report falls back to the thresholded matrix.
dir.create(empty <- file.path(tempdir(), "gistic_cont_empty"), showWarnings = FALSE)
stopifnot(is.null(load_gistic_continuous_at(empty, cnv)))
stopifnot(is.null(load_gistic_continuous_at(file.path(tempdir(), "no_such_dir_xyz"), cnv)))

cat("load_gistic_continuous: ALL PASS\n")
