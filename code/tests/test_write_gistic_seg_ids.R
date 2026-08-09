# write_gistic_seg(ids=) drives the per-group GISTIC runs. A silent mismatch here
# (e.g. ids that match nothing) would produce an EMPTY seg file and a GISTIC run on
# zero samples, so no-match must return NULL rather than an empty file.
source(file.path("code", "load_wes_results.R"))

dir.create(sd <- file.path(tempdir(), "segids"), showWarnings = FALSE)
for (s in c("S1", "S2", "S3")) {
  write.table(
    data.frame(Chromosome = "chr1", Start = 1, End = 1000,
               Num_Markers = 50, Segment_Mean = 0.5),
    file.path(sd, paste0(s, ".seg")), sep = "\t", quote = FALSE, row.names = FALSE)
}

cnv <- list(seg = list(dir = basename(sd), id_strip = "-1TAD104|_tumor_only",
                       chrom_col = "Chromosome", start_col = "Start", end_col = "End",
                       value_col = "Segment_Mean", num_col = "Num_Markers",
                       value_is_log2 = TRUE,
                       gistic_seg_out = file.path(tempdir(), "all.seg")))

out_all <- write_gistic_seg(cnv, out = file.path(tempdir(), "all.seg"))
stopifnot(!is.null(out_all), file.exists(out_all))
stopifnot(length(unique(read.delim(out_all)$Sample)) == 3)

out_sub <- write_gistic_seg(cnv, out = file.path(tempdir(), "sub.seg"),
                            ids = c("S1", "S3"))
stopifnot(!is.null(out_sub))
stopifnot(identical(sort(unique(read.delim(out_sub)$Sample)), c("S1", "S3")))

stopifnot(is.null(write_gistic_seg(cnv, out = file.path(tempdir(), "none.seg"),
                                   ids = "NOSUCHSAMPLE")))

cat("write_gistic_seg(ids=): OK\n")
