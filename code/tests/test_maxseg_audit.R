# GISTIC drops samples over -maxseg SILENTLY. Hypersegmentation correlates with
# aneuploidy, so the cut preferentially removes aneuploidy-high tumours — biasing
# exactly the axis under study. The audit must name every dropped sample and test
# whether exclusion is associated with group.
source(file.path("code", "attend_classes.R"))
source(file.path("code", "attend_scna.R"))

counts <- data.frame(ID = paste0("S", 1:6),
                     n_seg = c(100, 200, 50000, 300, 47000, 400))
groups <- data.frame(ID = paste0("S", 1:6),
                     scna_group = factor(c("MMRd-high", "MMRd-high", "MMRd-high",
                                           "MMRp-low", "MMRp-low", "MMRp-low"),
                                         levels = attend_scna$group_levels))

a <- maxseg_audit(counts, groups, maxseg = 46000)

stopifnot(identical(sort(a$excluded$ID), c("S3", "S5")))
stopifnot(a$n_before == 6, a$n_after == 4)
stopifnot(all(attend_scna$group_levels %in% rownames(a$tab)))  # empty levels kept
stopifnot(is.numeric(a$fisher_p))

# Nothing excluded -> empty frame, NA p, counts unchanged.
b <- maxseg_audit(data.frame(ID = "S1", n_seg = 10),
                  data.frame(ID = "S1",
                             scna_group = factor("MMRd-high",
                                                 levels = attend_scna$group_levels)),
                  maxseg = 46000)
stopifnot(nrow(b$excluded) == 0, is.na(b$fisher_p), b$n_after == 1)

cat("maxseg_audit: OK\n")
