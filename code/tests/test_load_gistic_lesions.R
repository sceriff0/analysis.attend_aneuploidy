# Synthetic all_lesions.conf_99.txt. GISTIC writes each peak TWICE — thresholded
# rows and "- CN values" rows. Keeping both would double every peak and inflate
# every frequency by 2x, so the parser must drop the CN-values rows.
source(file.path("code", "load_wes_results.R"))

dir.create(td <- file.path(tempdir(), "gistic_lesions"), showWarnings = FALSE)

les <- data.frame(
  `Unique Name` = c("Amplification Peak 1", "Deletion Peak 1",
                    "Amplification Peak 1 - CN values", "Deletion Peak 1 - CN values"),
  Descriptor = c("8q24.21", "10q23.31", "8q24.21", "10q23.31"),
  `Wide Peak Limits` = c("chr8:127735434-128753674(probes 100:200)",
                         "chr10:87863113-87971930(probes 300:400)",
                         "chr8:127735434-128753674(probes 100:200)",
                         "chr10:87863113-87971930(probes 300:400)"),
  `Peak Limits`   = rep("na", 4),
  `Region Limits` = rep("na", 4),
  `q values`      = c(1e-8, 1e-5, 1e-8, 1e-5),
  `Residual q values after removing segments shared with higher peaks` = rep(1e-4, 4),
  `Broad or Focal` = rep("Focal", 4),
  `Amplitude Threshold` = rep("0: t<0.3; 1: 0.3<t<0.9; 2: t>0.9", 4),
  `S1_tumor_only` = c(2, 0, 1.4, 0.0),
  S2              = c(0, 1, 0.0, -0.5),
  check.names = FALSE
)
write.table(les, file.path(td, "all_lesions.conf_99.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cnv <- list(gistic = list(dir = basename(td), all_lesions_glob = "*all_lesions*.txt"),
            seg = list(id_strip = "-1TAD104|_tumor_only"))

x <- load_gistic_lesions_at(td, cnv)
stopifnot(!is.null(x), all(c("peaks", "mat") %in% names(x)))

# CN-values rows dropped: 2 peaks, not 4.
stopifnot(nrow(x$peaks) == 2)
stopifnot(identical(x$peaks$direction, c("amp", "del")))
stopifnot(identical(x$peaks$chrom, c("8", "10")))
stopifnot(x$peaks$wide_start[1] == 127735434, x$peaks$wide_end[1] == 128753674)
stopifnot(x$peaks$q_value[1] == 1e-8)

# Matrix: samples x peaks, id_strip applied to sample names.
stopifnot(identical(sort(rownames(x$mat)), c("S1", "S2")))
stopifnot(ncol(x$mat) == 2)
stopifnot(x$mat["S1", x$peaks$peak_id[1]] == 2)
stopifnot(x$mat["S2", x$peaks$peak_id[2]] == 1)
stopifnot(is.integer(x$mat[1, 1]) || is.numeric(x$mat[1, 1]))

# Absent folder is knit-safe.
stopifnot(is.null(load_gistic_lesions_at(file.path(tempdir(), "nope"), cnv)))

cat("load_gistic_lesions_at: OK\n")
