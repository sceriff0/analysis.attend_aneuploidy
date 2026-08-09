# =============================================================================
# load_clinical.R  —  CLINICAL / IMAGING DATA LOADERS
#
# Sourced by reports 1, 3, 4, 5 and 7. Each function returns a data frame with
# the ID column the pipeline expects (see make_id_cfg() in attend_harmonise.R):
#   load_imaging_data()         -> ID (= pid) + image_id      [data/attend_data.xlsx]
#   load_clinical_data()        -> ID (= pid)                 [data/attend_data.xlsx]
#   load_gianlu_clinical_data() -> ID (= pid) + TUMOR_BARCODE [data/clinical_gianlu_data.csv]
# load_gianlu_clinical_data() is what builds the barcode->pid crosswalk.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(readxl)
})

# Shared workbook — both Excel-backed loaders read different sheets of this file.
attend_xlsx <- function() here("data", "attend_data.xlsx")

load_imaging_data <- function() {
  imaging_data <- read_excel(attend_xlsx(), sheet = "Imaging (2)")

  imaging_data |>
    rename(
      ID                     = `Unique Subject Identifier`,
      patient_id             = `Codice caso studio IEO`,
      image_id               = `file names`,
      MSI_Status             = `Evaluation of Microsatellite status, decode`,
      treatment              = `Treatment Arm`,
      Neoplastic_cellularity = `Neoplastic cellularity (%)`,
      Immune_infiltrate      = `Immune infiltrating cells (%)`,
      ARID1A                 = `ARID1A expression`,
      MMR_Status             = `Evaluation of MMR status`,
      Tumor_cells_PDL1       = `Percentage of tumor cells (TC) expressing PD-L1`,
      PFS_Months             = `PFS time`,
      PFS_Event              = `PFS event (PD/Death/no event)`
    ) |>
    rename_with(~ "Immune_infiltrate_PD1",  matches("expressing PD-1")) |>
    rename_with(~ "Immune_infiltrate_PDL1", matches("infiltrating immune cells")) |>
    select(
      ID, patient_id, image_id, MSI_Status, treatment,
      Neoplastic_cellularity, Immune_infiltrate, ARID1A, MMR_Status,
      Immune_infiltrate_PD1, Tumor_cells_PDL1, Immune_infiltrate_PDL1,
      PFS_Months, PFS_Event
    )
}

load_clinical_data <- function() {
  clinical_data <- read_excel(attend_xlsx(), sheet = "ClinData")

  clinical_data |>
    rename(
      ID                     = `Unique Subject Identifier`,
      IEO_ID                 = `Codice caso studio IEO`,
      MSI_STATUS             = `Evaluation of Microsatellite status, decode`,
      TREATMENT              = `Treatment Arm`,
      RACE                   = `Race`,
      MALIGNANCY             = `Malignancy`,
      HISTO_TYPE             = `Histological type`,
      HISTO_GRADE            = `Histological grade`,
      STATUS                 = `Status of disease`,
      BIOPSY_SITE            = `Type of sample provided for the analyses...12`,   # or ...13
      NEOPLASTIC_CELLULARITY = `Neoplastic cellularity (%)`,
      IMMUNE_INFILTRATE      = `Immune infiltrating cells (%)`,
      ARID1A                 = `ARID1A expression`,
      CLONAL_LOSS_AREA       = `If Clonal loss, percentage of involved tumor area (%)`,
      MLH1                   = `Mlh1 expression`,
      PMS2                   = `PMS2 expression`,
      MSH6                   = `Msh6 expression`,
      MSH2                   = `Msh2 expression`,
      MMR_STATUS             = `Evaluation of MMR status`,
      FIGO                   = `F.I.G.O. Stage`,
      PFS_EVENT              = `PFS event (PD/Death/no event)`,
      PFS_MONTHS             = `PFS time`
    ) |>
    rename_with(~ "IMMUNE_PD1",  matches("expressing PD-1")) |>
    rename_with(~ "IMMUNE_PDL1", matches("infiltrating immune cells")) |>
    select(
      ID, IEO_ID, MSI_STATUS, TREATMENT, RACE, MALIGNANCY, HISTO_TYPE,
      HISTO_GRADE, STATUS, BIOPSY_SITE, NEOPLASTIC_CELLULARITY, IMMUNE_INFILTRATE,
      ARID1A, CLONAL_LOSS_AREA, MLH1, PMS2, MSH6, MSH2, MMR_STATUS,
      IMMUNE_PD1, IMMUNE_PDL1, FIGO, PFS_EVENT, PFS_MONTHS
    )
}

load_gianlu_clinical_data <- function() {
  file_path <- here("data", "clinical_gianlu_data.csv")
  clinical_data <- read.csv(file_path, sep = "\t")

  clinical_data |>
    select(ID = Unique_Subject_Identifier, TUMOR_BARCODE = Tumor_Sample_Barcode, SEQ_CODE = Seq_name,
           MSI_STATUS = Evaluation.of.Microsatellite.status..decode,
           TREATMENT = Treatment.Arm, RACE = Race,  MALIGNANCY = Malignancy, HISTO_TYPE = Histological.type,
           HISTO_GRADE = Histological.grade, STATUS = Status.of.disease,
           NEOPLASTIC_CELLULARITY = Neoplastic.cellularity...., IMMUNE_INFILTRATE = Immune.infiltrating.cells....,
           ARID1A = ARID1A.expression, CLONAL_LOSS_AREA = If.Clonal.loss..percentage.of.involved.tumor.area....,
           MLH1 = Mlh1.expression, PMS2 = PMS2.expression, MSH6 = Msh6.expression, MSH2 = Msh2.expression,
           MMR_STATUS = Evaluation.of.MMR.status,  FIGO = F.I.G.O..Stage,
           PFS_EVENT = PFS.event..PD.Death.no.event., PFS_MONTHS = PFS_time)
}
