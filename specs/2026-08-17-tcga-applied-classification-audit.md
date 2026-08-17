# Audit: the ATTEND classification applied to TCGA data

**Date:** 2026-08-17
**Status:** two of three reports ported (10, 11); the third (below) audited, not ported.

## What was missing, and where it was

Three reports existed in the **sibling `ATTEND/` repository's git history** and never
reached `attend_aneuploidy`. They were deleted from `ATTEND` itself during the
consolidation to eight reports; the last commit that carries all three is
`70bc111` (`:construction: WIP baseline — reports 09i/09j/16/17, cross-cohort 12,
TCGA-2013 pipeline`).

| Deleted file (in `ATTEND@70bc111`) | Lines | Fate here |
|---|---|---|
| `analysis/12_mmrd_aneuploidy_crosscohort.Rmd` | 637 | **Ported** → `analysis/11-mmrd-crosscohort.Rmd` |
| `analysis/17_published_peak_clustering.Rmd` | 326 | **Ported** → `analysis/10-published-peak-clustering.Rmd` |
| `analysis/16_tcga2013_fig1a.Rmd` | 935 | **Not ported.** Audited below. |

Recovering a copy:

```bash
git -C ../ATTEND show 70bc111:analysis/16_tcga2013_fig1a.Rmd > /tmp/16_tcga2013_fig1a.Rmd
```

### This explains the orphaned-helper list

`CLAUDE.md` records four tested helpers with no caller. Three of the four are orphaned for
exactly one reason — the reports that used to call them are the three above:

| Helper | Called by |
|---|---|
| `mmrd_aneuploidy_split()` (`attend_classes.R` §15) | old report 12 → now report 11 |
| `load_tcga_2013_ascets()` (`load_wes_results.R`) | old report 12 → now report 11 |
| `load_tcga_published_peaks()` (`load_wes_results.R`) | old report 17 → now report 10 |
| `segments_to_bin_matrix()` (`attend_classes.R`) | old report 16 → **still orphaned** |

So `segments_to_bin_matrix()` is the one entry on that list which is still genuinely
without a caller after this pass, and porting report 16 is what would re-wire it.

## What report 16 did

**Title:** *"TCGA UCEC 2013 — Figure 1a reproduced, with aneuploidy."*

It ran **this pipeline's own clustering machinery on TCGA's own segments**, where the right
answer is published, and measured whether the pipeline recovers it. Two stated purposes:

1. **A methods control.** Report 09 runs `segments_to_bin_matrix()` → `cluster_arm_matrix()`
   → `fig1a_heatmap()` on ATTEND and asks which copy-number cluster is the serous-like one.
   Report 16 runs the identical pipeline on TCGA's segments, where TCGA ships its own
   assignment (`CNA_CLUSTER_K4`), and scores agreement with an **adjusted Rand index**. If
   the pipeline cannot recover TCGA's partition on TCGA's data, report 09's ATTEND clusters
   need a second look before they are interpreted.
2. **To place aneuploidy on that figure.** The 2013 paper has no per-sample aneuploidy
   score, so the relationship between its four SCNA clusters and aneuploidy was never drawn.
   Every heatmap in report 16 carries an aneuploidy annotation bar.

It drew **four partitions of the same cohort as four heatmaps** sharing one annotation
stack, so they can be read against each other:

1. **Published** — TCGA's own `CNA_CLUSTER_K4`; the reference answer.
2. **Genome-bin, Euclidean + Ward.D2** — the naive feature space: 1 Mb log2 ratios
   genome-wide, significant or not.
3. and 4. further feature-space / metric variants, swept and scored against (1).

## Overlap with what now exists

Report 10 (ported) already covers the *external-feature-set* half of this question: it
clusters TCGA on TCGA's **published** peaks and scores the result against `CNA_CLUSTER_K4`
with the same adjusted Rand index, using the pipeline's configured `cnv_dist`/`cnv_method`.

What report 10 does **not** cover, and only report 16 did:

- **The genome-bin feature space** (`segments_to_bin_matrix()`, 1 Mb bins, significant or
  not) as a control against peak-restricted features. This is the comparison that says
  whether restricting to significant peaks helps or hurts.
- **A distance/linkage sweep.** Report 10 deliberately uses the one configured pair so that
  any difference from report 09 is the feature set. Report 16 varied the metric.
- **Aneuploidy annotation on the TCGA Fig-1a heatmap**, i.e. where aneuploidy sits relative
  to TCGA's own four clusters.
- **A same-cohort GISTIC comparison** — report 16 also ran our GISTIC parameters on TCGA
  (`code/run_gistic_tcga.sh`) and compared. Report 10 flags this as an open question rather
  than answering it, because this repository has no such TCGA GISTIC run.

## What porting report 16 would take

1. **Rename and renumber** to `analysis/12-tcga2013-fig1a.Rmd`. Reports are numbered
   consecutively in reading order and 10/11 are now taken; 12 appends without renumbering.
2. **Scope contract.** Add the four-line Question / Cohort / Reads / Out-of-scope
   blockquote, and delete every other blockquote — `test_rmd_style.R` rule [4] permits only
   the contract. Report 16 carries several long "Reading it" blockquotes.
3. **Derivation notes.** One to three sentences of prose directly above each of its visible
   chunks, tracing to the raw data (`data/tcga_2013/data_cna_hg19.seg`), not to the nearest
   intermediate. Rule [2].
4. **Strip interpretation.** Report 16 ends in takeaway-style prose; conclusions belong in
   the manuscript. Statistical caveats move *up* into the derivation note of the chunk they
   qualify — the pattern used in reports 10 and 11.
5. **Figure system.** Replace `theme_minimal()`/`theme_bw()` with `attend_theme()`, remove
   every literal hex and base-R colour name (`"red"`, `"steelblue"`, `"turquoise3"`,
   `"firebrick"` all appear), and route any boxplot through `attend_box()`.
   `test_plot_style.R` rules [1]–[4] cover all of this.
6. **Line cap.** 935 lines against a `MAX_LINES` of 800. Either split the metric sweep out
   or add a named `OVER_CAP` entry with a reason. Preference: the four-heatmap comparison is
   one question and should stay together, so an `OVER_CAP` entry is the honest option —
   though trimming the prose that duplicates report 10 should get it under the cap on its
   own.
7. **`get_master()`.** Report 16 predates `code/build_master.R`; replace any
   `read_intermediate("attend_master_joined")` with `get_master()` so it inherits the
   no-build-order guarantee.
8. **Index entry** in `analysis/index.Rmd`, under **External**.

No missing library code. Every helper report 16 needs already exists here and is tested:
`segments_to_bin_matrix()`, `cluster_arm_matrix()`, `adjusted_rand_index()`,
`cnv_high_cluster_burden()`, `fig1a_heatmap()`, `feature_positions()`,
`load_tcga_ucec_2013()`, `load_tcga_2013_seg()`, `load_tcga_2013_ascets()`.

Data is present too: `data/tcga_2013/` holds `data_clinical_patient.txt`,
`data_clinical_sample.txt`, `data_cna_hg19.seg` and `pancan_aneuploidy.tsv`. The
shared-ASCETS panels of report 11 additionally need `bash code/run_ascets_tcga.sh` to have
been run; the report skips them and says so when it has not.
