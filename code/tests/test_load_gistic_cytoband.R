# Synthetic all_thresholded.by_genes.txt: col1 gene, col2 Locus ID, col3 Cytoband, then samples.
# NOTE: load_wes_results.R does `library(tidyverse)` at top — this test therefore requires
# tidyverse to be installed. In the bootstrap/dev env where tidyverse is absent, verify instead
# with: Rscript -e 'invisible(parse("code/load_wes_results.R")); cat("PARSE OK\n")' and defer
# this runtime check to the restored renv (see task-4-report.md).
source(file.path("code", "load_wes_results.R"))
dir.create(td <- file.path(tempdir(), "gistic_cb"), showWarnings = FALSE)
thr <- data.frame(`Gene Symbol` = c("MYC","TP53"), `Locus ID` = c(1,2), Cytoband = c("8q24.21","17p13.1"),
                  S1 = c(2, -2), S2 = c(1, 0), check.names = FALSE)
write.table(thr, file.path(td, "all_thresholded.by_genes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

cnv <- list(gistic = list(dir = basename(td),
                          all_thresholded_glob = "*all_thresholded.by_genes.txt",
                          amp_genes_glob = "*amp_genes*.txt", del_genes_glob = "*del_genes*.txt"),
            seg = list(id_strip = "-1TAD104|_tumor_only"))
x <- tryCatch(load_gistic_thresholded_at(td, cnv), error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
stopifnot(!is.null(attr(x, "feature_pos")))
fp <- attr(x, "feature_pos")
stopifnot(all(c("feature","cytoband") %in% names(fp)))
stopifnot("8q24.21" %in% fp$cytoband)
cat("load_gistic cytoband attr: OK\n")
