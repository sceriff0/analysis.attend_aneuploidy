# Base-R unit test. attend_harmonise.R does library(tidyverse) at its top; shim it
# so the pure-base-R recovery logic runs on a machine without tidyverse (delegates
# to real library() when a package IS installed — cluster-safe).
local({
  real_library <- base::library
  assign("library", function(package, ...) {
    pkg <- as.character(substitute(package))
    if (length(find.package(pkg, quiet = TRUE)) == 0) return(invisible())
    real_library(pkg, character.only = TRUE, ...)
  }, envir = globalenv())
})
source(file.path("code", "attend_harmonise.R"))

cfg <- list(col = "ID", space = "barcode")
cw  <- data.frame(barcode = c("ABC123", "XYZ789"), pid = c("P1", "P2"),
                  stringsAsFactors = FALSE)
img <- data.frame(image = "IMG1", pid = "P1", stringsAsFactors = FALSE)
df  <- data.frame(ID = c("ABC123", "abc123", " XYZ789 ", "NOPE"),
                  stringsAsFactors = FALSE)

## norm_id: upper + trim
stopifnot(identical(norm_id(c(" aB c ", "X")), c("AB C", "X")))

## recover = FALSE — exact only (today's behaviour, unchanged)
pv0 <- make_pid_vector(cw, img, recover = FALSE)
stopifnot(identical(pv0(df, cfg), c("P1", NA_character_, NA_character_, NA_character_)))

## recover = TRUE — case/whitespace near-misses promoted; genuine miss stays NA
pvr <- make_pid_vector(cw, img, recover = TRUE)
stopifnot(identical(pvr(df, cfg), c("P1", "P1", "P2", NA_character_)))

## collision guard — two barcodes normalise to the same key but map to DIFFERENT pids
## => that norm-key is ambiguous and must be REFUSED (stay NA), never guessed.
cwx <- data.frame(barcode = c("dup", "DUP"), pid = c("P1", "P2"), stringsAsFactors = FALSE)
pvx <- make_pid_vector(cwx, img, recover = TRUE)
stopifnot(is.na(pvx(data.frame(ID = "Dup", stringsAsFactors = FALSE), cfg)))

## same norm-key but SAME pid (two barcode formats for one patient) => safe, promoted
cws <- data.frame(barcode = c("s1", "S1"), pid = c("P9", "P9"), stringsAsFactors = FALSE)
pvs <- make_pid_vector(cws, img, recover = TRUE)
stopifnot(identical(pvs(data.frame(ID = "s1 ", stringsAsFactors = FALSE), cfg), "P9"))

## pid-space passthrough is unaffected by recover
stopifnot(identical(
  make_pid_vector(cw, img, recover = TRUE)(data.frame(ID = c("P1","P2")),
                                           list(col = "ID", space = "pid")),
  c("P1", "P2")))

cat("test_id_recovery: ALL PASS\n")
