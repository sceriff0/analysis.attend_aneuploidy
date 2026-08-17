# Where the paper-faithful TCGA Fig. 1a heatmap lives

**Date:** 2026-08-17
**Status:** diagnosed and fixed. Report 09 now displays the genome-wide continuous matrix
(commit `a68533f`); `analysis/10-tcga2013-fig1a.Rmd` carries the same figure on TCGA's own
data, which re-wires `segments_to_bin_matrix()` — the last helper on CLAUDE.md's orphan list.

## The complaint

Report 09's heatmap is **not the heatmap in Kandoth et al., *Nature* 2013, Fig. 1a**.

The paper's figure is: *"Tumours were hierarchically clustered into four groups based on
SCNAs. The heat map shows SCNAs in each tumour (horizontal axis) plotted by **chromosomal
location** (vertical axis)"* — a genome-wide, **continuous** copy-number landscape.

Report 09 plots something else, and it is one line that makes the difference
(`analysis/09-tcga-classification.Rmd:119`):

```r
cnv_mat <- tryCatch(load_gistic_thresholded(), error = function(e) NULL)  # ID x peak-gene {-2..+2}
```

and then `fig1a_heatmap(..., value_type = "thresholded")`. So the rows are **GISTIC
significant-peak genes** and the values are **discretised** {−2,−1,0,+1,+2}. Placing peak
genes by cytoband keeps the axis *ordered* by genome position but it is not the *continuous*
chromosomal axis of the paper — most of the genome is simply absent, and the colour scale is
five steps rather than a gradient.

**`fig1a_heatmap()` is not the problem.** It is byte-identical to `ATTEND` HEAD, and since it
was written (`ATTEND@3979409`) it has only gained an NA-cluster guard and a centralised
top-annotation builder. It already supports the faithful mode — `value_type = "continuous"` —
and `feature_positions()` already parses `chr1:0`-style genome-bin column names into
chromosome + coordinate. The machinery is intact; report 09 just does not feed it the right
matrix.

## When it was right, and when it was lost

Chronology from `ATTEND`'s history:

| Date | Commit | Event |
|---|---|---|
| 2026-07-12 | `7996157` | genome-binned Fig-1a heatmap first built from the DRAGEN `.seg` |
| **2026-07-12** | **`aa7a40d`** | **the regression.** *"GISTIC thresholded-peak clustering is the sole CN clustering in report 09"* — it *"removed the arm-level and bin-level clustering paths"*. From here report 09's heatmap is the peak matrix. |
| 2026-07-13 | `3979409` | `fig1a_heatmap()` extracted — *"faithful TCGA Fig-1a via ComplexHeatmap (chr-location rows, cluster-split tumours, red/blue SCNA)"* |
| **2026-07-13** | **`a77b0f9`** | **report `09c_tcga_classification_binned_k4.Rmd` — "CNV genome-bin cascade, k=4 (faithful TCGA Fig-1a + PFS KM)". This is the right plot.** |
| 2026-07-13 | `21f8ec1` | report `09d` — the same at k = 2 |
| 2026-07-13 | `9ab23f8`, `7d0e7dd` | reports `09e`/`09f` — arm-level "Fig-1a *variant*" (39 arms, explicitly not faithful) |
| 2026-07-28 | `70bc111` | last commit holding the deleted 09c/09d/16/17 |

Recover the faithful report:

```bash
git -C ../ATTEND show a77b0f9:analysis/09c_tcga_classification_binned_k4.Rmd > /tmp/09c.Rmd
```

## What 09c did that report 09 does not

```r
seg_long <- load_seg_data()                                  # data/seg/, DRAGEN .seg
arm_bnds <- load_arm_boundaries()                            # data/arm_boundaries_hg38.tsv
cnv_mat  <- segments_to_bin_matrix(seg_long, arm_bnds)       # sample x 1 Mb genome bin, CONTINUOUS
cl       <- cluster_arm_matrix(cnv_mat, k = k_cnv, method = "complete", dist_method = "euclidean")
fig1a_heatmap(cl$matrix, feat_pos, cl$cluster, ann_col = ann, value_type = "continuous", ...)
```

Its own summary of the difference, which is the argument for restoring it:

- **Continuous, not thresholded** — raw log2 seg-means, matching Fig. 1a's continuous colour
  scale rather than GISTIC's {−2..+2} calls.
- **Genome-wide, not peak-restricted** — every 1 Mb bin on every chromosome contributes, so
  broad SCNAs invisible to a peak-only feature space are shown, *"exactly as Fig. 1a shows
  the whole chromosomal landscape"*.
- **No GISTIC dependency** — built straight from the `.seg`, no external tool run, so it does
  not wait on `runbooks/gistic.md`.

## Everything needed is already in this repo

| Piece | Where | State |
|---|---|---|
| `load_seg_data()` | `code/load_wes_results.R` | present |
| `load_arm_boundaries()` | `code/attend_classes.R` | present |
| `segments_to_bin_matrix()` | `code/attend_classes.R` | present, **orphaned** — 09c was its only caller |
| `fig1a_heatmap(value_type = "continuous")` | `code/attend_plots.R` | present, unchanged since `3979409` |
| `feature_positions()` bin parsing | `code/attend_plots.R` | present (`^chr[0-9XY]+:`) |
| `data/arm_boundaries_hg38.tsv` | `data/` | present |
| `data/seg/` (DRAGEN `.seg`) | `data/seg/` | **empty in this checkout** — lives on HPC |

So this is a wiring job, not a rebuild. `segments_to_bin_matrix()` being on CLAUDE.md's
orphaned list is the same symptom seen from the library side.

## Two decisions before restoring

1. **Where it goes.** Report 09 currently declares GISTIC peaks *"the only clustering
   input"*, and `00-methods.Rmd` justifies that choice. The faithful heatmap can be an
   additional panel in 09 (peaks keep the cascade), can replace 09's clustering input
   (undoing `aa7a40d`), or can come back as its own report.
2. **The metric.** 09c clustered with **complete linkage + Euclidean** (R `hclust`
   defaults). The pipeline now mandates **1 − Pearson + Ward.D2** as its only
   distance/linkage pair (`cnv_dist`/`cnv_method` in `attend_classes.R`, justified in
   `00-methods.Rmd`). Reproducing 09c's figure exactly means using complete+Euclidean and
   documenting the exception; using the mandated pair means the same feature space but a
   different dendrogram from the figure that was liked.

## Also in ATTEND history, not restored

Deleted at the same consolidation, listed here so they are not lost again:

| File (`ATTEND@70bc111`) | Lines | What it was |
|---|---|---|
| `analysis/16_tcga2013_fig1a.Rmd` | 935 | The pipeline run on **TCGA's own** segments, scored against published `CNA_CLUSTER_K4` with an adjusted Rand index; four partitions as four heatmaps, aneuploidy annotation added. **Superseded 2026-08-17 by `analysis/10-tcga2013-fig1a.Rmd`**, which answers the same question with the two-matrix design (published peaks clustered, 1 Mb genome bins displayed) in 460 lines instead of 935 — the metric sweep is not reproduced, since the pipeline now mandates one distance/linkage pair |
| `analysis/17_published_peak_clustering.Rmd` | 326 | Both cohorts in TCGA's published 79-peak space |
| `analysis/12_mmrd_aneuploidy_crosscohort.Rmd` | 637 | ATTEND metastatic vs TCGA primary MMRd aneuploidy split |
| `analysis/09d_tcga_classification_binned_k2.Rmd` | — | genome-bin faithful Fig-1a at k = 2 |
| `analysis/09e`/`09f` | — | arm-level Fig-1a *variant*, k = 4 / k = 2 |
| `analysis/09g`/`09h`/`09i`/`09j` | — | GISTIC k = 2–8 sweeps, Pearson/Ward and Euclidean/Ward, directional |

17 and 12 were ported on 2026-08-17 and then withdrawn the same day — not needed for now.
