# Implementation plan — recurrent-CNA audit fixes

**Branch:** `cna-audit-fixes`  **Base:** `b0188e0`
**Companion spec:** `specs/2026-08-11-report-reorganization-design.md`

Fixes two confirmed denominator bugs in the recurrent copy-number figures, brings both
figures onto field convention, adds a memo-sorted oncoplot, and records the audit in the
design spec.

Both bugs are reproduced empirically (see Task 1). They are wrong in any report layout, so
this plan is independent of the 8-to-15 reorganization and lands before it.

---

## Global Constraints

**GC1 — Target toolchain is R 4.3.2 / Bioconductor 3.18.** Add no new package
dependencies. Everything here uses packages the reports already load: `dplyr`, `tidyr`,
`tibble`, `stringr`, `purrr`, `ggplot2`.

**GC2 — Knit-safe.** Every helper keeps its `NULL` in -> `NULL` out contract. Optional
dependencies stay behind `requireNamespace()`. A missing data folder must never error.

**GC3 — Config, not literals.** New thresholds go in `attend_cnv` in
`code/attend_classes.R`. Never hard-code a path, threshold, or column name at a call site.

**GC4 — Backwards compatible.** `recurrent_arm_calls()` must keep returning a list with
`long`, `freq`, `top_arms`, `wide_top`; `segment_pileup()` must keep returning a tibble
with `chrom`, `bin_start`, `gpos`, `gain`, `loss` plus the `chrom_offsets` and `n_samples`
attributes. Add fields; do not remove or rename existing ones. Existing callers
(`analysis/05_oncoplots_recurrent_cna.Rmd`, the oncoplot annotation path) must keep working.

**GC5 — Tests are pure-transform and run offline.** No real data, no `data.table`, no
`tidyverse` meta-package, no network. Build fixtures inline with `tibble::tibble()`.
Follow the style of `code/tests/test_scna_stats.R`: `source()` the file under test, then
`stopifnot()` on named properties.

**GC6 — Test command.** The local machine has R 4.6 and the tidyverse *member* packages
but not the meta-package; a shim library supplies it:

```bash
export SHIM=/private/tmp/claude-501/-Users-valer-Desktop-Github-attend-aneuploidy/cbe80688-99a4-45ac-aa2e-884f6312ecb5/scratchpad/rlib
cd /Users/valer/Desktop/Github/attend_aneuploidy-cna-fixes
R_LIBS_USER="$SHIM" Rscript code/tests/<test-file>.R && echo PASS
```

A test that prints nothing and exits 0 has passed. Report the exact command and its output.

**GC7 — Do not reformat untouched code.** Diffs stay minimal and reviewable.

**GC8 — Do not run `janitor::clean_names()`** on any prefixed/keyed table — it breaks
every `attend_cols` lookup.

---

## Task 1 — Fix the two denominator bugs (P0)

**Files:** `code/attend_classes.R`, new `code/tests/test_recurrent_cna_denominators.R`

Both bugs are confirmed. This reproduction currently FAILS on both counts and must PASS
when you are done:

```r
source("code/attend_classes.R")

# BUG 1 — 4 samples, 1 arm: 2 AMP, 1 NEUTRAL, 1 LOWCOV (arm not evaluable)
calls <- tibble::tibble(ID = paste0("S", 1:4), `1q` = c("AMP","AMP","NEUTRAL","LOWCOV"))
recurrent_arm_calls(calls)$freq$gain
#   currently 0.50  (LOWCOV counted in the denominator)
#   correct   0.667 (2 of 3 evaluable arms)

# BUG 2 — 4 samples profiled, only 2 carry an altered segment
segs <- tibble::tibble(ID = c("S1","S2"), Chromosome = "chr1",
                       Start = c(0,0), End = c(5e5,5e5), Type = "DUP", logRatio = 1.0)
arms <- tibble::tibble(chrom = "chr1", arm = c("1p","1q"),
                       start = c(0,1e6), end = c(1e6,2e6))
segment_pileup(segs, arms, bin = 1e6, min_abs = 0)$gain[1]
#   currently 1.0  (denominator = samples WITH alterations)
#   correct   0.5  (denominator = samples PROFILED), given n_samples = 4
```

### 1a — `recurrent_arm_calls()` (currently at `code/attend_classes.R:548`)

ASCETS emits `LOWCOV` for any arm below the breadth-of-coverage threshold. ASCETS's own
aneuploidy score — and this repo's `load_wes_results.R:485` — divide by
`sum(call != "LOWCOV")`. The recurrence barplot must use the same denominator.

Change `norm_call()` from a 3-way `case_when` with a `TRUE ~ "neutral"` catch-all to an
explicit 4-way classification:

| Tokens (upper-cased, trimmed) | Maps to |
|---|---|
| `AMP`, `GAIN`, `1`, `+1`, `2` | `"gain"` |
| `DEL`, `LOSS`, `-1`, `-2` | `"loss"` |
| `NEUTRAL`, `NEUT`, `0` | `"neutral"` |
| `LOWCOV`, `NC`, `NA`, `""`, and `NA` itself | `NA_character_` |
| anything else | `NA_character_`, and `warning()` once listing the distinct unrecognised tokens |

The catch-all must NOT silently become `"neutral"` — that is the bug.

In the `freq` summarise, keep `na.rm = TRUE` and add:

- `n_evaluable = sum(!is.na(call))` — the denominator, reported per arm
- `n_lowcov    = sum(is.na(call))`

Add a `min_evaluable_frac` argument (default `0.5`). An arm evaluable in fewer than that
fraction of samples is still reported in `freq`, but is **excluded from `top_arms`** — a
p-arm evaluable in 2 of 40 samples must not become an oncoplot annotation track on the
strength of a 2-sample denominator. Rank `top_arms` by `altered` among eligible arms only.

`freq` keeps its existing columns (`arm`, `gain`, `loss`, `altered`) and stays sorted by
`desc(altered)`.

### 1b — `segment_pileup()` (currently at `code/attend_classes.R:586`, denominator at `:616`)

`n_samples <- dplyr::n_distinct(segs$ID)` runs *after* the `min_abs` and direction
filters, so a sample whose genome is quiet contributes no rows and leaves the denominator.

Add an `n_samples` argument, default `NULL`. When `NULL`, compute it from `cnv_long`
**before any filtering**: `dplyr::n_distinct(cnv_long[[id_col]])`. When supplied, use the
given value (so a caller can pass the full cohort size, covering samples whose annotation
file failed to parse). Never derive it from the filtered `segs`.

Keep `attr(wide, "n_samples")` reporting the denominator actually used.

### Tests — `code/tests/test_recurrent_cna_denominators.R`

Cover, with `stopifnot()` and a comment naming the property:

1. LOWCOV excluded from the arm denominator (the 2/3 case above).
2. `NC`, empty string, and `NA` all behave as no-call.
3. An unrecognised token maps to `NA` and raises a warning.
4. `n_evaluable` / `n_lowcov` are correct.
5. `min_evaluable_frac` keeps a poorly-covered arm out of `top_arms` but in `freq`.
6. Pileup denominator counts profiled samples, not altered ones (the 0.5 case above).
7. An explicit `n_samples` argument overrides the derived value.
8. `NULL` in -> `NULL` out still holds for both functions (GC2).
9. Return shapes still carry every field named in GC4.

---

## Task 2 — Genome ordering and amplitude tiers (P1)

**Files:** `code/attend_classes.R`, new `code/tests/test_cna_plot_inputs.R`

Field convention (Beroukhim 2010; Zack 2013; GenVisR `cnFreq`) is that recurrence is
displayed **along the genome**, and that low-level gain/loss is separated from high-level
amplification/homozygous deletion.

### 2a — `order_arms_genomic(arms)`

New exported helper in `code/attend_classes.R`. Takes a character vector of arm names
(`"1p"`, `"1q"`, ..., `"22q"`, `"Xp"`, `"Xq"`, optionally `chr`-prefixed) and returns them
ordered `1p, 1q, 2p, 2q, ... 22p, 22q, Xp, Xq, Yp, Yq`. Unparseable names sort last, in
their original relative order. Case-insensitive; tolerates a `chr` prefix.

Return the ordered character vector. Do not mutate `recurrent_arm_calls()`'s existing sort
— the report chooses the order (GC4).

### 2b — amplitude tiers in `segment_pileup()`

Add a `high_abs` argument, defaulting to a new `attend_cnv$pileup$high_abs_logratio`
config entry set to `1.0` (|log2 ratio| >= 1, roughly a full copy gained or lost — the
conventional low-level/high-level boundary; document the choice in a comment).

Return four new columns alongside the existing `gain` / `loss`:

- `gain_low`, `gain_high` — fraction of profiled samples whose gain segment overlapping
  that bin is below / at-or-above `high_abs`
- `loss_low`, `loss_high` — the same for losses

`gain` and `loss` keep their current meaning: the fraction with **any** gain / loss
(GC4). Where a sample has both a low and a high segment in the same bin and direction,
count it in the **high** tier only, and once — no sample may be double-counted within a
(bin, direction).

Segments whose `logRatio` is `NA` count toward the low tier.

### Tests — `code/tests/test_cna_plot_inputs.R`

1. `order_arms_genomic()` puts `10p` after `9q` and before `11p` (string sort would not).
2. It tolerates `chr` prefixes and mixed case, and sorts X/Y last.
3. Unparseable names sort last without erroring.
4. A |logRatio| = 2 segment lands in `gain_high`, not `gain_low`.
5. `gain == gain_low + gain_high` when no sample has both tiers in one bin.
6. A sample with both a low and a high gain segment in one bin counts once, in `gain_high`
   — so `gain` still equals `1/n_samples`, not `2/n_samples`.
7. `NA` logRatio counts as low tier.
8. `NULL` in -> `NULL` out (GC2); GC4 fields all still present.

---

## Task 3 — Recurrent-CNA figures on convention (P1/P2)

**File:** `analysis/05_oncoplots_recurrent_cna.Rmd` only. No changes under `code/`.

The report currently draws both figures with `theme_minimal(base_size = 11)`
(lines ~274, ~319, ~359) instead of the pipeline's `attend_theme()`, sorts arms by
frequency, and shows a single binary gain/loss tier.

If `source(here("code", "attend_plots.R"))` is not already among the setup chunk's
sources, add it — `attend_theme()` lives there.

### 3a — arm-level barplot (chunk `recurrent-cnv`)

- Order arms by `order_arms_genomic()`, not by `desc(altered)`. Keep the diverging
  gain-right / loss-left layout and `coord_flip()`.
- Use `attend_theme()`.
- Print each arm's `n_evaluable` — either as an axis-label suffix (`1q (n=38)`) or a
  companion column in the printed table. A reader must be able to see that an arm's
  frequency rests on a reduced denominator.
- Update the prose above the chunk: state that the denominator is evaluable arms
  (LOWCOV excluded) and that arms are in genome order.

### 3b — genome-wide pileup (chunk `segment-pileup`)

- **Two amplitude panels**, low-level above and high-level below, sharing an x axis
  (`facet_grid(tier ~ ., scales = "free_y")` is acceptable). Gains up and red, losses
  down and blue, as now.
- **Centromere markers:** a light dashed vertical line at each chromosome's p/q boundary,
  read from the arm-boundary table already passed in (`load_arm_boundaries()`), not
  hard-coded.
- **Drop chrY.** The cohort is uterine and therefore all female; a chrY axis segment is
  an empty tail. Filter it out at the plotting step, with a one-line comment saying why.
  Keep chrX.
- Use `attend_theme()`.
- Update the prose: name the denominator (samples profiled), the amplitude threshold, and
  state plainly that **these are descriptive frequencies with no significance testing** —
  GISTIC q-values are the field's significance standard and are not applied here. One
  sentence; do not write an essay (the spec's prose rules apply).

Both chunks keep their existing `eval` / `is.null()` guards so the report still knits
without `data/ascets/` or `data/cnv_annotations/` (GC2).

---

## Task 4 — Memo-sorted oncoplot as the primary panel

**File:** `analysis/05_oncoplots_recurrent_cna.Rmd` only.

The report draws its oncoplots with `sortByAnnotation = TRUE` (chunks `oncoplot` at
~line 170 and `oncoplot-panel` at ~line 189). That is a deliberate, documented choice, and
it is the right figure for showing enrichment — but it overrides the memo sort that is the
cBioPortal/OncoPrint default and the layout readers expect as the primary mutation
landscape.

Add a **new chunk before the existing `oncoplot` chunk**, drawing the same top-20 gene
oncoplot with the default memo sort (`sortByAnnotation = FALSE`), the same
`clinicalFeatures = covar`, and `removeNonMutated = FALSE` (whole-cohort denominators —
do not change this).

Give it a heading that makes the pair legible, e.g.
`## Oncoplot — mutation landscape (memo-sorted)` for the new one, and retitle the existing
one `## Oncoplot — sorted by annotation (enrichment view)`.

Two or three sentences of prose above the new chunk: memo sort orders samples so
co-occurring and mutually exclusive alterations form the characteristic staircase; the
annotation-sorted panel below answers a different question. Do not restate what an
oncoplot is.

Keep the existing `eval = MAFTOOLS_OK` gate. Do not touch the `mutfreq-groups` chunk.

---

## Task 5 — Record the audit in the design spec

**File:** `specs/2026-08-11-report-reorganization-design.md` only.

Append a new section `## 11. Recurrent-CNA audit (2026-08-11)` after the existing
`## 10. Open items`, so the document's numbering stays monotonic. Do not reflow or reword
any existing section.

The section records:

1. **The two denominator bugs**, each with: the file and line as of `b0188e0`, the
   observed vs correct value from Task 1's reproduction, and the internal contradiction
   (`load_wes_results.R:485` and `run_ascets_tcga.sh:12` already use the evaluable-arm
   denominator; `05:170`'s `removeNonMutated = FALSE` already insists on whole-cohort
   denominators for mutations).
2. **The convention gaps** — genome ordering, amplitude tiers, centromere marks, the
   `theme_minimal` inconsistency, chrY — as a table with the source for each convention.
3. **The absence of significance testing**, and that GISTIC output already exists in
   `data/gistic/`.
4. **The oncoplot finding** — memo sort vs annotation sort, and that both are now drawn.
5. A **P0/P1/P2 table** mapping each change to the file it lands in and to its destination
   report under the Config B split (`34_recurrent_cna`, `33_mutation_landscape`).
6. Sources, as markdown links:
   - GISTIC 2.0 — https://broadinstitute.github.io/gistic2/
   - GenePattern GISTIC_2.0 — https://www.genepattern.org/modules/docs/GISTIC_2.0/6.3/
   - ASCETS — https://github.com/beroukhim-lab/ascets
   - GenVisR cnFreq — https://genviz.org/module-03-genvisr/0003/05/01/cnFreq_GenVisR/
   - Zack et al. 2013 — https://www.nature.com/articles/ng.2760
   - TCGA SCNA analysis, Front. Oncol. 2021 — https://www.frontiersin.org/articles/10.3389/fonc.2021.700568/full
   - cBioPortal OncoPrinter — https://www.cbioportal.org/oncoprinter
   - ComplexHeatmap OncoPrint — https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html

Also update §8's phase table so the audit fixes appear as work that precedes the split, and
§9 ("Out of scope") so it no longer implies the CNA figures are untouched.

Prose follows the spec's own §5 rules: at most two bold spans per paragraph, numbers as
facts not adjectives, no takeaway section.

---

## Task 6 — Make the corrected pileup denominator reach production (P0 follow-up)

**Added after Task 1.** Task 1 implemented §1b exactly as written, and verification showed
the specification was incomplete: the corrected denominator is opt-in and nothing opts in.

**Files:** `code/load_wes_results.R`, `code/attend_classes.R`,
new `code/tests/test_pileup_profiled_denominator.R`

### The gap

Task 1 made `segment_pileup()` count distinct ids in `cnv_long` before filtering. That
fixes the case where a sample has segments that all fail the `min_abs` cut. It cannot fix
the case where a sample contributes **no rows at all** — and that is the production case:

- `load_cnv_data()` (`code/load_wes_results.R:174`) reads one `*.cnv.annotated.tsv` per
  sample. ClassifyCNV annotates *called* CNVs, so a copy-number-quiet tumour yields a file
  with zero data rows, and `mutate(ID = id)` on a zero-row tibble yields zero rows. The
  sample disappears before `segment_pileup()` ever sees it.
- The only call site, `analysis/05_oncoplots_recurrent_cna.Rmd:334`, is
  `segment_pileup(load_cnv_data(), load_arm_boundaries())` — no `n_samples` argument.

Verified against the Task 1 code:

| Fixture | `n_samples` used | `gain` reported | Correct |
|---|---|---|---|
| 4 samples, 2 filtered out by `min_abs` | 4 | 0.5 | 0.5 — fixed |
| 4 samples, 2 contribute no rows | 2 | 1.0 | 0.5 — still wrong |
| same, with explicit `n_samples = 4` | 4 | 0.5 | 0.5 — correct |

### 6a — `load_cnv_data()` reports how many samples it profiled

In `code/load_wes_results.R`, `load_cnv_data()` already computes `files`, one path per
sample. Attach the count to the returned tibble:

```r
attr(out, "n_profiled") <- length(files)
```

Set it on every return path that reads a folder, including the path that returns an empty
tibble because the folder has no matching files (`n_profiled = 0`). The early return for a
missing folder keeps returning a bare `tibble()` with no attribute — absence of the folder
is not the same as zero profiled samples.

Preserve the attribute through the final `relocate()` — dplyr verbs drop non-standard
attributes, so set it **after** the last dplyr call, immediately before returning.

### 6b — `segment_pileup()` prefers the profiled count

In `code/attend_classes.R`, when the `n_samples` argument is `NULL`, resolve the
denominator in this order:

1. `attr(cnv_long, "n_profiled")`, when present, non-`NA`, and at least the number of
   distinct ids actually present. (The floor guard means a stale or wrong attribute can
   never *shrink* the denominator below the samples in hand.)
2. otherwise `dplyr::n_distinct(cnv_long[[id_col]])`, computed before any filtering, as
   Task 1 implemented.

An explicit `n_samples` argument still wins over both. Record which source was used in a
new attribute `attr(wide, "n_samples_source")`, one of `"argument"`, `"profiled"`,
`"observed"`, so the report can state its denominator honestly.

### 6c — Tests, `code/tests/test_pileup_profiled_denominator.R`

`load_cnv_data()` needs `data.table`, which is **not installed on this machine**, so 6a is
not locally testable. Test 6b only, by setting the attribute on an inline fixture:

1. A fixture carrying `attr(x, "n_profiled") <- 4` with only 2 ids present gives
   `gain = 0.5`, `attr(., "n_samples") == 4`, `n_samples_source == "profiled"`.
2. No attribute present falls back to the observed distinct-id count;
   `n_samples_source == "observed"`.
3. An explicit `n_samples` argument overrides a present attribute;
   `n_samples_source == "argument"`.
4. The floor guard: an attribute of 1 with 2 ids present does not shrink the denominator —
   the observed count wins.
5. `NA` attribute is ignored, falling back to observed.
6. `NULL` in -> `NULL` out still holds (GC2).

State in the report that 6a is untested locally and why.

### 6d — Report prose

Do not edit `analysis/05_oncoplots_recurrent_cna.Rmd` in this task — Task 3 owns that file
and must state the denominator and its source in the pileup prose. Note the requirement in
your report so Task 3 picks it up.
