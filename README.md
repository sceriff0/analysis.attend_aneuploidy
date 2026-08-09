# attend_aneuploidy

Aneuploidy-focused multi-omic analysis of the ATTEND endometrial cancer cohort, built
with [workflowr](https://workflowr.github.io/workflowr/) and
[renv](https://rstudio.github.io/renv/).

- **R 4.3.2**, **Bioconductor 3.18** — pinned versions; `renv` is not yet initialised
  here (see *Status* below).

## Quick start

```r
renv::restore()              # once renv.lock exists
workflowr::wflow_build()     # knit every report into docs/
workflowr::wflow_build("analysis/02_survival_stratified.Rmd")   # knit one
```

## The pipeline

Reports run **in order**. Report 1 builds the per-patient master table
(`output/clean_data/attend_master_joined`); reports 2–8 read it back and analyse it.

| # | Report | Cohort |
|---|---|---|
| 1 | Data integration — the join, and who it leaves behind | all |
| 2 | Survival, stratified by molecular axis | all |
| 3 | Aneuploidy vs mutations, one gene at a time | **MMR-deficient only** |
| 4 | IHC/imaging concordance + cell-type composition | Pt 1 all imaged; **Pt 2 MMRd only** |
| 5 | Mutation landscape + recurrent copy number + aneuploidy×MMR contrast | all |
| 6 | Covariate distributions by MMR status | all |
| 7 | TCGA integrated classification (Pearson/Ward, directional serous call) | all |
| 8 | Is the TMB distribution split? (and does TCGA agree) | all |

Two conventions hold everywhere: **TMB is always exactly two definitions** — *normal*
(`tmb__TMB_SCORE`) and *ancestry corrected* (`tmb__TMB_nassar`) — and **highlighted
points** always mean the same thing (gold = has IHC/imaging data; magenta = sample
21S188, *polipo*).

See [`analysis/index.Rmd`](analysis/index.Rmd) for the full index and
[`CLAUDE.md`](CLAUDE.md) for the architecture.

## Status

Implemented: the workflowr layout, all eight reports, and the `code/*.R` helper layer.
The `code/load_*.R` loaders are real — they read the HPC exports.

Outstanding:

1. `renv` is not initialised — run `renv::init(bioconductor = "3.18")` **under R 4.3.2**.
2. The HPC inputs (`data/attend_barcodes.csv`, `data/seg/`, `data/variant_annotations/`,
   `data/ascets/`, `data/gistic/`) are not in this checkout, so reports knit to a
   guarded skeleton until they are synced.
3. Several `attend_cols` entries in `code/attend_classes.R` are placeholders and must be
   set to the real column names before results are valid.

## What is here but not in git

Deliberately untracked (see `.gitignore`): `papers/` (17 MB of reference PDFs),
`research/`, `data/tcga_2013/` (the TCGA UCEC reference cohort report 8 needs), and
`output/inventory/` (62 MB, regenerable via `code/inventory_data.R`).
