# Layer 3 (peak-set structure) and §4.2 (LOO stability). The asymmetric rule the
# whole report rests on: PRESENCE in the MMRd-high set is informative; ABSENCE is
# not, because absence is indistinguishable from low power at n=9. These helpers
# must therefore report retention FRACTIONS and support counts, never a bare
# present/absent verdict.
source(file.path("code", "attend_classes.R"))
source(file.path("code", "attend_scna.R"))

u <- data.frame(
  source   = c("pooled", "pooled", "mmrd_high", "mmrp_high"),
  union_id = c("U001_8_amp", "U002_10_del", "U001_8_amp", "U003_19_amp"),
  stringsAsFactors = FALSE
)
st <- peak_set_structure(u)
stopifnot(setequal(names(st$sets), c("pooled", "mmrd_high", "mmrp_high")))
stopifnot(st$jaccard["pooled", "mmrd_high"] == 0.5)   # {U1,U2} vs {U1}
stopifnot(st$jaccard["mmrd_high", "mmrp_high"] == 0)
stopifnot("U003_19_amp" %in% st$private$union_id[st$private$source == "mmrp_high"])
stopifnot(!"U001_8_amp" %in% st$private$union_id)     # shared -> not private

# --- profile correlation ----------------------------------------------------
set.seed(1)
m <- rbind(hi1 = c(1, 1, 0, 0), hi2 = c(1, 1, 0, 0),
           lo1 = c(0, 0, 1, 1), lo2 = c(0, 0, 1, 1),
           q   = c(1, 1, 0, 0))
g <- factor(c("MMRp-high", "MMRp-high", "MMRp-low", "MMRp-low", "MMRd-high"),
            levels = attend_scna$group_levels)
pc <- profile_correlation(m, g, target = "MMRd-high",
                          refs = c("MMRp-high", "MMRp-low"))
stopifnot(nrow(pc) == 1, pc$ID == "q")
stopifnot(pc$best_match == "MMRp-high")
stopifnot(pc[["MMRp-high"]] > pc[["MMRp-low"]])

cat("peak_set_structure + profile_correlation: OK\n")
