# add_scna_group(): crosses is_mmrd() with aneuploidy_class into a 4-level factor.
# NA MMR status must propagate to NA (not silently become MMRp) — that bug would
# reclassify every unknown-IHC tumour into the comparison group.
source(file.path("code", "attend_classes.R"))

df <- data.frame(
  MMR_class       = c("Deficient", "Deficient", "Intact", "Intact", NA),
  aneuploidy_class = c("aneuploidy-high", "aneuploidy-low",
                       "aneuploidy-high", "aneuploidy-low", "aneuploidy-high"),
  stringsAsFactors = FALSE
)

out <- add_scna_group(df)

stopifnot(all(c("mmr_group", "scna_group") %in% names(out)))
stopifnot(is.factor(out$scna_group))
stopifnot(identical(levels(out$scna_group),
                    c("MMRp-low", "MMRp-high", "MMRd-low", "MMRd-high")))
stopifnot(identical(as.character(out$scna_group),
                    c("MMRd-high", "MMRd-low", "MMRp-high", "MMRp-low", NA)))
stopifnot(identical(as.character(out$mmr_group),
                    c("MMRd", "MMRd", "MMRp", "MMRp", NA)))

# Panel config: 12 directional loci, 8 amp + 4 del.
stopifnot(nrow(attend_scna$panel) == 12)
stopifnot(sum(attend_scna$panel$direction == "amp") == 8)
stopifnot(sum(attend_scna$panel$direction == "del") == 4)
stopifnot(all(c("locus", "cytoband", "direction") %in% names(attend_scna$panel)))
stopifnot(attend_scna$maxseg == 46000)

cat("add_scna_group + attend_scna config: OK\n")
