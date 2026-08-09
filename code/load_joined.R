library(here)
library(fs)

clinical_code <- here("code", "load_clinical.R")
ihc_code <- here("code", "load_phenotypes.R")
wes_code <- here("code", "load_wes_results.R")

conversion_file <- here("data", "attend_barcodes.csv")
joined_file <- here("data", "joined_data.csv")

source(clinical_code)
source(ihc_code)
source(wes_code)


if (!file_exists(joined_file)) {
  clinical_data <- load_clinical_data()
  imaging_data <- load_imaging_data()
  ihc_data <- load_ihc_data()
  aneu_data <- load_aneuploidy_scores()
  maf_data <- load_maf_data()
  conversion_data <- read.csv(conversion_file)
} else {
  joined_data <- read.csv(joined_file)
}

dim(clinical_data)
dim(ihc_data)
dim(aneu_data)
dim(maf_data)
dim(conversion_data)

data <- read.csv(here("data", "lets.csv"), sep = "\t")
dim(data)
