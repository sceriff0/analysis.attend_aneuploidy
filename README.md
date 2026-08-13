# attend_aneuploidy

Aneuploidy-focused multi-omic analysis of the ATTEND endometrial cancer cohort, built
with [workflowr](https://workflowr.github.io/workflowr/) and
[renv](https://rstudio.github.io/renv/).

- **R 4.3.2**, **Bioconductor 3.18** — pinned versions; `renv` is not yet initialised
  here (see *Status* below).

## Quick start

```r
renv::restore()              # once renv.lock exists
workflowr::wflow_build()     # knit every report into docs/ — no build order to remember
```

There is **no build order**. The join lives in `code/build_master.R`; every report calls
`get_master()`, which reads the master table when present and builds it when not. Any
report can be knitted first. To rebuild the master by hand:

```bash
Rscript code/build_master.R
```

```bash
Rscript code/tests/test_rmd_style.R   # report structure (base R, no deps)
Rscript code/tests/test_rmd_parse.R   # every chunk parses
```

## The pipeline

Nine documents, numbered in reading order `00` to `08`. Each opens with a four-line scope
contract — Question / Cohort / Reads / Out of scope.

| # | Report | Cohort |
|---|---|---|
| 00 | Methods — every justification and shared definition, once | — |
| 01 | Data integration — the join, and what it cost | all |
| 02 | FlowPath vs the pathologist | all imaged |
| 03 | Covariate distributions by MMR status | all |
| 04 | Survival — every KM and Cox in the pipeline | all |
| 05 | Aneuploidy in MMR-deficient tumours — mutations, then composition | **MMR-deficient only** |
| 06 | Mutation landscape and copy number — oncoplots, recurrent CNA, subgroup contrast | all |
| 07 | Tumour mutational burden — bimodality, by aneuploidy x MMR, TCGA replication | all |
| 08 | TCGA integrated classification (Pearson/Ward, directional serous call) | all |

Shared definitions live once in `analysis/00_methods.Rmd` and are **linked** from the
reports, not rendered into each one. Operational recipes live in `runbooks/`. Reports carry
no takeaway sections — conclusions belong in the manuscript.

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
