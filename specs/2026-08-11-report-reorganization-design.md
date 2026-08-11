# Report reorganization — design

**Date:** 2026-08-11
**Status:** approved, not yet implemented
**Scope:** `analysis/*.Rmd` structure and prose. No change to `code/*.R` behaviour.

---

## 1. Problem

### 1.1 Reports mix unrelated analyses

Five of eight reports are structurally more than one document. Each announces this itself
with a `# Part N` header.

| Report | Lines | Contains |
|---|---|---|
| `01_data_integration` | 785 | the join + master write (Part 1); the unmapped/coverage audit (Part 2); value sanity (Part 3) |
| `04_ihc_imaging_spatial` | 806 | FlowPath-vs-pathologist concordance on all imaged patients (Part 1); cell-type composition vs aneuploidy in MMRd only (Part 2) |
| `05_oncoplots_recurrent_cna` | 657 | oncoplots + subgroup mutation frequency (Part 1); recurrent CNA (Part 2); a pre-specified permutation-tested subgroup contrast (Part 3) |
| `07_tcga_classification` | 1020 | the TCGA cascade (46–618); a Nassar TMB methods essay (619–921); a GISTIC SLURM runbook (922–1005) |
| `08_tmb_distribution` | 840 | bimodality battery (179–412); TCGA UCEC replication (413–620); TMB by aneuploidy within MMR strata (621–810) |

Topic-occurrence counts (grep hits per report) show the smear:

```
                                    TMB ANEU MMR SURV TP53 GIST ASCE  IHC TCGA ONCO
01_data_integration                  37   11    0    1    1    0    0   11    0    0
02_survival_stratified               22   31   42   26    0    0    1    0    0    0
03_aneuploidy_mutations_mmrd         21   40   52    0   12    0    1    1    0    0
04_ihc_imaging_spatial                9  107   69   16   15    0    0   77    0    0
05_oncoplots_recurrent_cna           15   18   49    3    3   40    7    1   15   41
06_ecdf_covariates                   13    1   26    0    0    0    0    0    0    0
07_tcga_classification              121   58   23    5   37   46   13   10   71    1
08_tmb_distribution                 167   73   90    1    0    0    1    0   71    0
```

### 1.2 Two confirmed consequences of the mixing

1. **Duplicated analysis.** `03` §5 ("Aneuploidy by TMB class", lines 288–343) and `08` §6
   ("TMB by aneuploidy status, within MMR strata", lines 621–810) are the same
   association drawn from opposite axes, on different cohorts. `03`'s version is
   MMRd-only and therefore silently drops the MMRp comparison; `08`'s stratifies by MMR
   and is strictly more informative.

2. **Scattered survival.** `02` is the survival report, yet `04` runs two further KM
   analyses (`survival-by-ihc` at 04:328, `survival-aneu-immune` at 04:774) because
   `immune_class` happens to be defined there. Analyses drift to wherever their exposure
   variable was computed.

### 1.3 Prose mixes three genres

| Genre | Example | Belongs |
|---|---|---|
| Derivation note — what is shown, from which raw file/columns, which interpretation-changing operations | `08:634–641` | inline, above the chunk |
| Justification — why Pearson/Ward, why the directional serous call, how Nassar works | `07:53–75`, `07:619–921` | `00_methods.Rmd` |
| Operations — how to run GISTIC under SLURM | `07:922–1005` | `runbooks/` |

`07` is 1020 lines because it contains all three. In-chunk comment density ranges from
2% (`01`, `04`) to 13% (`07`) — there is no house rule.

Specific density problems: `07:53–75` carries 11 bold spans in one blockquote; `05:443`
hardcodes "n ≈ 9" inside a *methodological* argument, so the justification for choosing a
composite score over per-peak tests will silently outlive the cohort size that motivates it.

### 1.4 The shared explainers are retyped

The "two TMB definitions" block appears in six of eight reports, and again in
`index.Rmd`. Highlight-group and cohort-restriction explanations repeat similarly.

---

## 2. Decisions

| # | Decision | Rationale |
|---|---|---|
| D1 | One `.Rmd` answers **one question on one cohort**. | If a report needs "Part 1 … Part 2" to describe itself, it is two reports. |
| D2 | **Banded numbering** `1x` build / `2x` QC / `3x` findings / `4x` external. | The number states the stage; gaps absorb future reports without renumbering. `ec_analysis/analysis/` already uses this scheme. |
| D3 | Shared explainers become **child fragments** in `analysis/_explainers/`, included via `child=`. | Written once; the text still renders into every report, so a single knitted HTML remains self-contained. Prose stays in Markdown rather than R strings. |
| D4 | **All survival consolidates into `30_survival`**, including the two KM chunks currently in `04`. | One question, one report. Cost: `30_survival` gains a dependency on the IHC join for `immune_class`. Accepted. |
| D5 | Justification prose moves to `00_methods.Rmd`; operations prose moves to `runbooks/`. | Both are read once, not per-knit. |
| D6 | **No takeaway/interpretation sections.** The tier is removed from the standard and the eight existing sections are deleted. | Reports show derivations and results; conclusions belong in the manuscript, and editorialised summaries drift as data changes. |
| D7 | `runbooks/` and `specs/` live at **repo root**, never under `docs/`. | workflowr owns `docs/` as its build output; `wflow_publish()` overwrites it. |

### D6 — one carve-out

Some content in the current `## What to take away` sections is a **statistical caveat**,
not a takeaway. Example, `06:163`:

> The KS test is uncorrected and the two TMB rows are the same patients twice.

Caveats of this kind are **not deleted** — they move up into the derivation note of the
chunk they qualify (Tier 2), where a reader meets them before the figure rather than after.
Only forward-looking interpretation is removed.

---

## 3. Target layout

```
analysis/
  00_methods.Rmd                     every justification, once. Knits the Nassar fit.
  _explainers/                       child fragments (not built by wflow_build)

  ── 1x  BUILD ────────────────────────────────────────────────────────
  10_data_integration.Rmd            join, collapse, TMB recompute. WRITES master.
  11_join_audit.Rmd                  unmapped records, coverage, UpSet, value sanity

  ── 2x  QC / VALIDATION ──────────────────────────────────────────────
  20_flowpath_concordance.Rmd        FlowPath vs pathologist        [all imaged]
  21_covariates_by_mmr.Rmd           ECDF + KS + boxplots by MMR    [all]

  ── 3x  FINDINGS ─────────────────────────────────────────────────────
  30_survival.Rmd                    all KM/Cox                     [all]
  31_aneuploidy_mutations.Rmd        aneuploidy vs gene status      [MMRd]
  32_cell_composition.Rmd            composition vs aneuploidy      [MMRd]
  33_mutation_landscape.Rmd          oncoplots + subgroup frequency [all]
  34_recurrent_cna.Rmd               arm + genome-wide CNA          [all]
  35_subgroup_cna_contrast.Rmd       MMRd-high vs MMRp-high         [contrast]
  36_tmb_distribution.Rmd            bimodality battery             [all]
  37_tmb_by_aneuploidy_and_mmr.Rmd   TMB boxplots, def x MMR group  [by MMR stratum]

  ── 4x  EXTERNAL ─────────────────────────────────────────────────────
  40_tcga_classification.Rmd         the cascade                    [all]
  41_tcga_tmb_replication.Rmd        TCGA UCEC replication          [TCGA]

runbooks/
  gistic.md                          07:922–964
  classification_options.md          07:965–1005
specs/
  2026-08-11-report-reorganization-design.md   (this file)
```

15 knitted analysis reports + `00_methods`, versus 8 today.

### 3.1 Scope contract

Every report opens with a fixed four-line block:

```markdown
> **Question.** Does cell-type composition differ by aneuploidy class?
> **Cohort.** MMR-deficient only (n = `r nrow(dat)`).
> **Reads.** `attend_master_joined` + the FlowPath annotation tree.
> **Out of scope.** Concordance with the pathologist -> [20](20_flowpath_concordance.html).
```

The `Out of scope` line is the enforcement mechanism. A contributor adding a TMB panel to
the composition report must either delete a claim or put the panel where it belongs.

---

## 4. Chunk destination map

Complete for content chunks. Of the 148 chunks, 133 carry content and each has exactly one
destination below. The remaining 15 are per-report `setup` / `session-info` boilerplate,
regenerated in each new report rather than moved. Nothing is dropped.

| Destination | Source chunks (report: labels, line range) |
|---|---|
| `10_data_integration` | 01: `load` … `master-write`, `pointblank-gate`, `save-each` (37–417) |
| `11_join_audit` | 01: `unmapped` … `msi-score-vs-status` (418–782) |
| `20_flowpath_concordance` | 04: `load`, `survival-join`, `ihc-metrics`, `concordance-cellularity`, `concordance-immune`, `concordance-pd1`, `concordance-aridia`, `concordance-pdl1`, `concordance-tcpdl1`, `ihc-patient` (104–327) |
| `21_covariates_by_mmr` | 06: all (`data`, `ecdf`, `ks`) + **new** aneuploidy-by-MMR boxplot (see §7.1) |
| `30_survival` | 02: all · **+ 04 `survival-by-ihc` (328), `survival-aneu-immune` (774)** |
| `31_aneuploidy_mutations` | 03: `restrict-mmrd` … `other-classes-box` (125–287) |
| `32_cell_composition` | 04: `data` (358) … `imaging-covars-cat` (358–773) |
| `33_mutation_landscape` | 05: `build-maf`, `maf-summary`, `oncoplot`, `oncoplot-panel`, `mutfreq-groups` (88–297) |
| `34_recurrent_cna` | 05: `recurrent-cnv`, `segment-pileup` (298–383) |
| `35_subgroup_cna_contrast` | 05: `scna-gate` … `scna-loo` (384–654) |
| `36_tmb_distribution` | 08: `bimod-helpers` … `synthesis-note` (127–434) |
| `37_tmb_by_aneuploidy_and_mmr` | 08: `aneu-tmb-build`, `aneu-tmb-box`, `aneu-tmb-test`, `aneu-tmb-note` (642–810) **+ 03 `tmb-definitions` (288–343), merged** |
| `40_tcga_classification` | 07: `cnv-load` … `cascade-panels-tmb` (148–618) |
| `41_tcga_tmb_replication` | 08: `tcga-load` … `tcga-mmr-contrast` (435–620) |
| `00_methods` | 07: `ancestry-tmb`, `ancestry-tmb-by-group`, `nassar-fit`, `nassar-effect` (722–921) + all justification prose extracted from 03/05/07 |
| `runbooks/gistic.md` | 07: 922–1005 (prose + bash, no R chunks) |
| *deleted* | the 8 `## What to take away` / `## Interpretation` / `## Notes` sections (141 lines), less the caveats promoted per D6 |

### 4.1 The 03 §5 merge

`03:288–343` (`tmb-definitions`) does not move — it **dissolves** into `37`. Its figure is
the MMRd facet of `37`'s `facet_grid(definition ~ mmr_group)`; `08`'s version already
draws both MMR groups. The merge removes a duplicated analysis rather than relocating it.

---

## 5. Prose standard

### 5.1 Three tiers

| Tier | Where | Length | Rule |
|---|---|---|---|
| 1. Scope contract | top of report, once | 4 lines | fixed shape (§3.1) |
| 2. Derivation note | above **every** visible-output chunk | 1–3 sentences | what is shown / where it comes from in the raw data / which operations change interpretation — **required** |
| 3. In-chunk comment | inside chunks | 1–2 lines | only non-obvious *code* decisions |

Tier 2 follows the `clean-rmd` rule already adopted by this project: trace to the *raw*
input (`data/ascets/`, `data/gistic/`), not to the nearest intermediate object.

**Banned inline:** justification (-> `00_methods.Rmd`), operations (-> `runbooks/`),
explainers repeated across reports (-> `_explainers/`), takeaway sections (D6).

### 5.2 Density rules

1. **At most 2 bold spans per paragraph.** `07:53–75` has 11.
2. **No blockquotes in `analysis/[1-4]*.Rmd`.** Every current one is justification.
   Blockquotes remain legal in `00_methods.Rmd` and in the scope contract.
3. **Every quoted number is inline R.** Not `n ≈ 9` but `` `r sum(grp == "MMRd-high")` ``.
4. **One clause per em-dash, at most one em-dash per sentence.**

### 5.3 Worked example

Before — `03:274–287`, 14 lines, half of it re-explaining TMB:

> Section 4's TMB panel uses the upstream score. Because ATTEND is sequenced tumour-only,
> residual germline variants inflate TMB in an ancestry-dependent way (Nassar et al.,
> *Cancer Cell* 2022), so the comparison is redrawn here under **both** definitions the
> pipeline carries: … Both are binarised at the same cutoff, so the only difference …

After — explainer becomes a child chunk; what remains is a derivation note:

    ```{r child='_explainers/tmb-two-definitions.Rmd'}
    ```

    Aneuploidy score (`aneu__aneuploidy_score` — ASCETS arm-level SCNA burden from
    `data/ascets/`) for TMB-high vs TMB-low patients, one facet per TMB definition.
    TMB is split at `r attend_thresholds$tmb` mut/Mb, so the facets differ only in who
    crossed the cut once the germline correction was applied. n = `r nrow(dat)`
    MMR-deficient patients; those without a TMB value are dropped.

    ```{r tmb-definitions, fig.height=5}

---

## 6. Explainer fragments

`analysis/_explainers/` is not built by `wflow_build()` (it globs `analysis/*.Rmd`,
non-recursively) and is ignored by the rmarkdown site generator because of the leading
underscore. Paths in `child=` resolve relative to `analysis/`.

| Fragment | Content | Used by |
|---|---|---|
| `tmb-two-definitions.Rmd` | the two reported definitions, the two banned intermediates, why tumour-only inflates TMB | 21, 30, 31, 32, 33, 36, 37 |
| `highlight-groups.Rmd` | `polipo` (21S188, magenta) and `cohort` (gold, multi-modal subset); first-match-wins ordering | 31, 32, 37 |
| `aneuploidy-classes.Rmd` | `attend_thresholds$aneuploidy` cut, high/low derivation, ASCETS provenance | 21, 31, 32, 34, 35, 37 |
| `cohort-mmrd.Rmd` | what "MMRd only" means, `is_mmrd()`, that MSI-high is a correlate not a criterion | 31, 32 |
| `master-provenance.Rmd` | what `attend_master_joined` is, written by `10`, read via `read_intermediate()` | all except 10 |

---

## 7. Enforcement

`code/tests/test_rmd_style.R`, alongside the existing suite:

- every `analysis/[0-4]*.Rmd` has a scope contract within its first 20 lines
- every chunk with visible output has at least one prose line directly above it
- no `analysis/*.Rmd` exceeds 400 lines
- no `^>` blockquote in `analysis/[1-4]*.Rmd`, except the four-line scope contract within
  the first 20 lines (§3.1)
- no `^## What to take away`, `^## Interpretation`, `^## Notes` anywhere
- the "two TMB definitions" table appears in exactly one place (`_explainers/`)

Parse-only; runs without `tidyverse`, so it passes on the local machine (which has R 4.6
and no tidyverse) as well as on the cluster.

### 7.1 New content

One genuine gap, found while mapping. The original
`attend_wes_analysis/analysis/00_get_aneuplodity_score_ascets_mod.Rmd:127–193` compares
**aneuploidy score itself between MMRd and MMRp** as an overlaid histogram plus ECDF with
the cutoff marked. In the current pipeline that comparison survives only as an ECDF in
`06` and a 2×2 contingency table at `05:432` — there is **no boxplot with a test**.
Given aneuploidy is this repo's subject, `21_covariates_by_mmr` gains that panel.

---

## 8. Migration order

Each phase leaves the site knittable.

| Phase | Work | Verify |
|---|---|---|
| audit | Recurrent-CNA audit (§11): three denominator bugs, genome ordering, amplitude tiers, centromere marks, chrY drop, memo-sorted oncoplot. Implemented on `cna-audit-fixes`, lands before phase 0. | `test_recurrent_cna_denominators.R`, `test_cna_plot_inputs.R`, `test_pileup_profiled_denominator.R`, `test_classifycnv_file_id.R` all pass |
| 0 | Create `_explainers/`, `00_methods.Rmd` (empty shell), `runbooks/`, `code/tests/test_rmd_style.R` (expected to fail). | test runs, reports still knit |
| 1 | Extract `07:619–1005` into `00_methods.Rmd` + `runbooks/`. No renumbering. | `07` drops to ~620 lines; knits |
| 2 | Write the five explainer fragments; replace the six in-report copies with `child=` chunks. | knit diff shows identical rendered text |
| 3 | Split the compounds: `01`->`10`,`11`; `04`->`20`,`32`; `05`->`33`,`34`,`35`; `08`->`36`,`37`,`41`; `07`->`40`. | all chunks accounted for against §4 |
| 4 | Renumber `02`->`30`, `03`->`31`, `06`->`21`. Update `_site.yml`, `index.Rmd`, `CLAUDE.md`. | no dead links |
| 5 | Merge `03` §5 into `37`; move the two KM chunks from `04` into `30` (D4). | `37` renders both MMR groups; `30` renders all KM |
| 6 | Prose pass: scope contracts, derivation notes, delete takeaway sections (promoting caveats per D6), density rules. Add the §7.1 boxplot. | style test green |

Phases 1 and 2 deliver most of the readability gain and require no renumbering, so they
can ship independently of 3–6.

---

## 9. Out of scope

- **Extracting analysis logic into `code/`.** `08:128–178` is a private bimodality library
  (skewness, kurtosis, Sarle's coefficient, mclust wrapper) living inside a report, and
  `04` Part 1 and `07`'s cascade have similar candidates. Worth doing, but after the
  boundaries settle. Tracked separately.
- `renv` initialization, the missing `data/attend_barcodes.csv`, and the placeholder
  `attend_cols` entries. Unrelated to this refactor.
- Any change to `code/*.R` behaviour by *this reorganization*. Only
  `code/tests/test_rmd_style.R` is added by it. The recurrent-CNA denominator and
  convention fixes described in §11 are a separate, prior effort on `cna-audit-fixes`
  and do touch `code/attend_classes.R` and `code/load_wes_results.R`.

## 10. Open items

- `04:328 survival-by-ihc` moves to `30_survival` per D4, which means `30` must load the
  IHC join. If that proves awkward at implementation time, the fallback is to keep it in
  `20` and note the exception in `30`'s `Out of scope` line.
- The 400-line cap in §7 is a first guess. `40_tcga_classification` may land above it even
  after the methods extraction; adjust the cap rather than force an unnatural split.

## 11. Recurrent-CNA audit (2026-08-11)

Before the split in §§3–8, report 05's two recurrent-CNA figures were audited against
field convention on branch `cna-audit-fixes` (base `b0188e0`). Three denominator bugs were
found and fixed, five convention gaps were closed, and a memo-sorted oncoplot was added.
No statistical test was added — that stays open. This section is the record.

### 11.1 Denominator bugs — one root, three instances

Each bug counted a *proxy* for "samples profiled" (rows that survived a filter, or files
found on disk) instead of the concept itself. Filtering and file-globbing are means of
finding profiled samples, not definitions of the term, and each proxy diverged from it in
its own report.

| # | Function | Location (`b0188e0`) | Proxy counted | Observed | Correct |
|---|---|---|---|---|---|
| 1 | `recurrent_arm_calls()` | `code/attend_classes.R:548`, catch-all at `:556` | rows not explicitly gain/loss/neutral | 0.500 (2 AMP, 1 NEUTRAL, 1 LOWCOV) | 0.667 (2 of 3 evaluable) |
| 2 | `segment_pileup()` | `code/attend_classes.R:586`, denominator at `:616` | rows surviving the magnitude/direction filters | 1.000 (4 profiled, 2 altered) | 0.500 |
| 3 | `load_cnv_data()` | `code/load_wes_results.R:174` | files found on disk | 3 files → `n_profiled = 3` for 2 distinct samples | 2 |

Bug 1: ASCETS emits `LOWCOV` for an arm below its breadth-of-coverage threshold. The
`case_when` catch-all `TRUE ~ "neutral"` swept `LOWCOV` into the denominator as evidence of
no alteration. `code/load_wes_results.R:485` and `code/run_ascets_tcga.sh:12` already
divided by `sum(call != "LOWCOV")` on the same file — the aneuploidy score and the
recurrence barplot used different denominators on identical input.

Bug 2: `n_samples <- dplyr::n_distinct(segs$ID)` ran after the magnitude and direction
filters, so a sample with no surviving altered segment vanished from the denominator
instead of counting as profiled-but-quiet. The oncoplot in the same report already passes
`removeNonMutated = FALSE` for exactly this reason — to keep whole-cohort denominators for
mutations.

Bug 3 surfaced while verifying the fix for bug 2, not in the original audit. ASCETS/
ClassifyCNV filenames get stripped to a sample id via `attend_cnv$classifycnv$id_strip`
(`-1TAD104|_tumor_only`), by design collapsing a resequenced sample's two filenames onto
one id — so `S1-1TAD104.cnv.annotated.tsv` and `S1_tumor_only.cnv.annotated.tsv` both
derive `S1`. Counting `length(files)` treats those as two profiled samples. This bug
inflates rather than deflates the denominator — the mirror of bug 2. The floor guard added
for bug 2 cannot catch it, since that guard only stops the attribute from being smaller
than the observed count, never larger.

All three are fixed. `recurrent_arm_calls()` gained an explicit four-way classification
(no catch-all), `n_evaluable`/`n_lowcov` columns, and a `min_evaluable_frac` gate that keeps
a poorly-covered arm out of `top_arms` while leaving it in `freq`. `segment_pileup()` now
resolves its denominator in order: an explicit `n_samples` argument, then
`attr(cnv_long, "n_profiled")` when present and no smaller than the observed count, then
the observed distinct-id count. It records which source won in
`attr(wide, "n_samples_source")`. `load_cnv_data()` attaches `n_profiled` as the count of
distinct ids derived by `.classifycnv_file_id()`, not `length(files)`.

### 11.2 Convention gaps closed

| Convention | Fix | Source |
|---|---|---|
| Arm display order | genome order (1p→22q, X, Y) via `order_arms_genomic()`; selection of the top 20 arms by alteration frequency is unchanged | [GenVisR `cnFreq`](https://genviz.org/module-03-genvisr/0003/05/01/cnFreq_GenVisR/) |
| Amplitude tiers | low/high split at `attend_cnv$pileup$high_abs_logratio = 1.0` (\|log2 ratio\| >= 1) | Beroukhim 2010; [Zack et al. 2013](https://www.nature.com/articles/ng.2760) |
| Centromere marks | dashed vertical line at each chromosome's p/q boundary, read from `load_arm_boundaries()` | [TCGA SCNA analysis, Front. Oncol. 2021](https://www.frontiersin.org/articles/10.3389/fonc.2021.700568/full) |
| chrY | dropped from the genome-wide pileup axis; the cohort is uterine and therefore all female; chrX is kept | cohort composition |
| Theme | `theme_minimal()` replaced by `attend_theme()` in both figures | pipeline house style |

### 11.3 No significance testing

Both figures remain descriptive frequencies; no test was added. GISTIC 2.0 q-values are
the field's standard for recurrence significance, and loaders for GISTIC output already
exist in this codebase (`code/attend_classes.R`, `code/load_wes_results.R`, consumed by
report 07's TCGA cascade from `data/gistic/`). Wiring GISTIC significance into the
recurrent-CNA figures themselves is **not done** — it stays outstanding.

Sources: [GISTIC 2.0](https://broadinstitute.github.io/gistic2/),
[GenePattern GISTIC_2.0](https://www.genepattern.org/modules/docs/GISTIC_2.0/6.3/),
[ASCETS](https://github.com/beroukhim-lab/ascets).

### 11.4 Oncoplot — memo sort restored

`sortByAnnotation = TRUE` on the existing `oncoplot` chunk is a deliberate, documented
choice for showing enrichment, but it overrode maftools' default memo sort — the
[cBioPortal OncoPrinter](https://www.cbioportal.org/oncoprinter) /
[ComplexHeatmap OncoPrint](https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html)
convention readers expect as the primary landscape. A new chunk, `oncoplot-memo`, draws the
same top-20 genes with `sortByAnnotation = FALSE`, the same `clinicalFeatures = covar`, and
`removeNonMutated = FALSE`, ahead of the existing chunk. Both panels now render: the
memo-sorted view first, the annotation-sorted view retitled as the enrichment view below it.

### 11.5 Verification limits

`code/load_wes_results.R` cannot execute in this development environment: `data.table`,
`fs`, and `here` are all unavailable. So bug 3's fix and bug 2's production wiring rest on
reading the code plus a test that pins `.classifycnv_file_id()` and `segment_pileup()`'s
resolution order, not a run of `load_cnv_data()` itself. The `.Rmd` changes (§11.2, §11.4)
were parse-checked only, via a manual chunk extractor and base `parse()`, since neither the
cohort data nor `knitr` is present locally. These carry weaker evidence than the helper
changes, which have passing unit tests; treat them accordingly, not as equally verified.

Recommended first check on the cluster: compare the on-disk count of
`*.cnv.annotated.tsv` files against `attr(load_cnv_data(), "n_profiled")`. A mismatch means
the id collision behind bug 3 is live in the real data.

### 11.6 Changes landed, by priority

Destinations are the target layout in §3 and the chunk destination map in §4.

| Priority | Change | File | Destination |
|---|---|---|---|
| P0 | `recurrent_arm_calls()`: LOWCOV excluded from the denominator, `n_evaluable`/`n_lowcov` added, `min_evaluable_frac` gates `top_arms` | `code/attend_classes.R` | `34_recurrent_cna` |
| P0 | `segment_pileup()`: denominator computed from `cnv_long` before filtering | `code/attend_classes.R` | `34_recurrent_cna` |
| P0 | `load_cnv_data()` attaches `n_profiled`; `segment_pileup()` prefers it over the observed count, floor-guarded, recording `n_samples_source` | `code/load_wes_results.R`, `code/attend_classes.R` | `34_recurrent_cna` |
| P0 | `n_profiled` derived from distinct `.classifycnv_file_id()` ids, not `length(files)` | `code/load_wes_results.R` | `34_recurrent_cna` |
| P1 | `order_arms_genomic()`; amplitude tiers (`gain_low`/`gain_high`/`loss_low`/`loss_high`) in `segment_pileup()` | `code/attend_classes.R` | `34_recurrent_cna` |
| P1 | Arm barplot: genome order, `n_evaluable` axis labels, `attend_theme()` | `analysis/05_oncoplots_recurrent_cna.Rmd` | `34_recurrent_cna` |
| P1/P2 | Genome-wide pileup: two amplitude panels, centromere marks, chrY dropped, `attend_theme()`, denominator and source stated in prose | `analysis/05_oncoplots_recurrent_cna.Rmd` | `34_recurrent_cna` |
| P2 | Memo-sorted oncoplot added as the primary panel; annotation-sorted panel retitled and kept | `analysis/05_oncoplots_recurrent_cna.Rmd` | `33_mutation_landscape` |
