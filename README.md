# attend_aneuploidy

Aneuploidy-focused multi-omic analysis of the ATTEND endometrial cancer cohort, built
with [workflowr](https://workflowr.github.io/workflowr/) and
[renv](https://rstudio.github.io/renv/).

- **R 4.3.2**, **Bioconductor 3.18** — pinned versions; `renv` is not yet initialised
  here (see *Status* below).

## Quick start

```r
renv::restore()                                    # once renv.lock exists
workflowr::wflow_build("analysis/10_data_integration.Rmd")  # FIRST: writes the master
workflowr::wflow_build()                           # then knit the rest into docs/
workflowr::wflow_build("analysis/30_survival.Rmd") # knit one

Rscript code/tests/test_rmd_style.R                # report structure (base R, no deps)
Rscript code/tests/test_rmd_parse.R                # every chunk parses
```

## The pipeline

Reports run **in order**, and the number states the stage: `1x` build, `2x` QC, `3x`
findings, `4x` external. Report 10 builds the per-patient master table
(`output/clean_data/attend_master_joined`); every later report reads it back.

Each report answers **one question on one cohort** and opens with a four-line scope
contract — Question / Cohort / Reads / Out of scope.

| # | Report | Cohort |
|---|---|---|
| 00 | Methods — every justification in the pipeline, once | — |
| 10 | Data integration — the join, the collapse, the TMB recompute | all |
| 11 | Join audit — unmapped records, coverage, value sanity | all |
| 20 | FlowPath vs the pathologist | all imaged |
| 21 | Covariate distributions by MMR status | all |
| 30 | Survival — every KM and Cox in the pipeline | all |
| 31 | Aneuploidy vs mutation status, one gene at a time | **MMR-deficient only** |
| 32 | Cell-type composition vs aneuploidy | **MMR-deficient only** |
| 33 | Mutation landscape — oncoplots + subgroup frequency | all |
| 34 | Recurrent copy number — arm level + genome-wide | all profiled |
| 35 | Aneuploidy x MMR copy-number contrast | MMRd-high vs MMRp-high |
| 36 | Is the TMB distribution split? | all |
| 37 | TMB by aneuploidy, once per MMR stratum | by MMR stratum |
| 40 | TCGA integrated classification (Pearson/Ward, directional serous call) | all |
| 41 | TCGA UCEC replication of the TMB split | TCGA |

Three reports write intermediates others read, so build order matters on a cold start:
**10** writes the master, and **36** writes `attend_tmb_long` / `attend_tmb_battery`,
which **37** and **41** read back. Alphabetical order satisfies this, except that `00`
sorts first and reads the master — so build `10` before the rest.

Shared explanations live once in `analysis/_explainers/` and are rendered into every
report that needs them via `child=`. Justification lives in `analysis/00_methods.Rmd`;
operational recipes live in `runbooks/`. Reports carry no takeaway sections — conclusions
belong in the manuscript.

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
