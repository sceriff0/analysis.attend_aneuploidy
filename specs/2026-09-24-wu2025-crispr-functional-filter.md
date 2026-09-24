# Porting Wu et al. (Immunity 2025)'s functional CNV filter onto ATTEND

**Status:** design, not implemented. No pipeline code changed.
**Date:** 2026-09-24
**Source paper skill:** `~/.claude/skills/wu-2025-immunity-paper` (built this session,
labelled **draft** — see §7).

---

## 0. One-paragraph summary

Wu et al. take a recurrence list of copy-number-altered genes and intersect it with a
catalogue of CRISPR-screen **resister** and **sensitizer** genes, direction-matched:
recurrently *deleted* ∩ *resisters*, recurrently *amplified* ∩ *sensitizers*. That
intersection is the functional filter that turns several thousand recurrently altered genes
into a finding about immune evasion. The gene catalogue is already converted and tracked in
this repo (`data/crispr_t_cell_screens_wu2025.tsv`, commit `58f13a1`) and **nothing consumes
it**. This spec proposes what should consume it, and — more importantly — records the three
places where a faithful transplant of their method would be *wrong* for ATTEND, and what to
do instead.

---

## 1. What Wu et al. actually did

Wu M, Yang S, Yang Z, …, Johnson DB, Liu S, Lo RS. "Genomic copy-number variants drive
apoptotic evasion underlying acquired resistance to immune checkpoint inhibitors."
*Immunity* 58:2864–2877 (2025). doi:10.1016/j.immuni.2025.10.001
`[FULL: complete main text and STAR Methods read page-by-page against rendered previews by
four independent reviewers over exclusive page ranges, plus the supplementary workbook mmc4
verified against its own OOXML; MISSING — Figures S1–S7 and Tables S1, S2, S4–S7, which were
never downloaded, so the membership of the 108/90 overlap sets and the per-patient recurrent
gene lists are UNVERIFIED]`

### 1.1 Their recurrence definition is a count, not a statistic

> "For the current clinical cohort, we defined recurrent DP-specific gene amplifications or
> deletions as ≥ three patients." — STAR Methods, printed e5

That is the entire definition. **The string "GISTIC" appears nowhere in the paper** — not in
the prose, and not in the Key Resources Table software list (checked independently by two
reviewers). There is no peak caller, no q-value, no amplitude threshold, no background model
and no significance test behind the word "recurrent".

- Segments: **union** of Sequenza 2.1.2 and VarScan2 2.4.3 calls on **WES** (238× mean).
- Purity/ploidy: Sequenza, default parameters.
- **No copy-number cutoff defining "amplification" or "deletion" is stated anywhere** for the
  clinical bulk cohort.
- **Reference build is never stated** for the clinical cohort (hg19 is stated only for the
  cell-line/model WES, WGS and scWGS). hg19 is an inference.
- "DP-specific" is a **set difference** against the patient's own pre-treatment baseline, with
  no CN delta, no log2 threshold and no clonality rule. (The quantitative ±0.09 ratio rule in
  the methods is **scWGS-only** and does not govern Figure 1.)

### 1.2 The gene catalogue

> "We defined sensitizer or resister genes as those hits with concordant evidence in ≥ two
> independent measurements." — STAR Methods, printed e5

23 published CRISPR-Cas9 screens of tumour-cell killing by T cells, CAR-T excluded. Per-hit
significance is delegated wholesale: *"according to the specific cutoffs defined by each
study."* Note "independent measurements" is **not** defined as "independent studies" — one
study with two arms appears to satisfy it.

### 1.3 The intersection, and what is not attached to it

Direction-matched: recurrently **deleted** ∩ **resisters** (loss of a resister promotes
resistance), recurrently **amplified** ∩ **sensitizers**.

**Figures 1D and 1E carry no test.** They are plain Venn counts: no p value, no background
universe, no multiple-testing correction. The hypergeometric test named in the Figure 1 legend
belongs to the Metascape REACTOME enrichment in panels 1B/1C — a different analysis on a
different gene set. Benjamini–Hochberg appears exactly once in the whole paper, inside
maftools' own `somaticInteractions`; MutSig2CV is run *deliberately* uncorrected.

Reported numbers: 519 resisters, 877 sensitizers, overlaps **108** and **90**; validation
cohort 117 and 277.

### 1.4 Three internal inconsistencies in the source, verified and unresolved

| Quantity | Conflict |
|---|---|
| Recurrence threshold | Results narrative **≥3** patients; Fig 1B/1C legend and Enrichment methods **≥4**; validation cohort **>2** |
| Deleted-gene count | Results says **5,035** in one sentence and **5,053** two sentences later |
| Cohort denominator | STAR Methods "patients, **n = 18**"; Figure 1 legend "of **17** patients" |

The 5,035/5,053 conflict is settled by reading the delivered Figure 1D crop directly: the Venn
is **4927 | 108 | 411**, and 4927 + 108 = **5,035**. Use 5,035; 5,053 is a typo.
The ≥3/≥4 and 17/18 conflicts are **not** resolvable from the material we have.

### 1.5 No code, restricted data

> "This paper does not report original codes." — Data and code availability, printed 2874

Bulk BAMs at EGA `EGAS50000001055`, scWGS at dbGaP `phs004102.v1`; **MAF files only from the
lead contact on request.** Anything built here is a reimplementation from prose, and the page
must say so.

---

## 2. What the repo already has

| Artifact | State |
|---|---|
| `code/fetch_crispr_screens.R` | **Implemented and correct.** Verified independently this session. |
| `data/crispr_t_cell_screens_wu2025.tsv` | **Tracked**, 5,789 rows. Verified. |
| `attend_crispr$min_measurements` | **Does not exist.** Named in the converter's header as the intended config knob. |
| Loader | **Does not exist.** |
| Any consumer | **None.** No report or helper references the file. |

### 2.1 The converter was audited against the workbook and passes

An independent OOXML re-read of `mmc4.xlsx` (zipfile + ElementTree, shared strings resolved by
index) compared **23,540 populated cells**: 0 coordinates present in only one, 0 value
mismatches, 141 merge ranges set-identical, 0 formula cells. Every count the converter claims
reproduces exactly:

| | Workbook (independent) | Repo TSV |
|---|---|---|
| resisters unfiltered / at ≥2 | 2,828 / **519** | 2,828 / **519** |
| sensitizers unfiltered / at ≥2 | 2,961 / **877** | 2,961 / **877** |
| DP-overlap YES | 108 / 90 | 108 / 90 |
| symbols in both lists at ≥2 | 105 | 105 |

The converter's positional read (`skip = 2`, explicit column indices 10/11/13 and 17/18/20) is
what makes it correct. The sheet is **four** blocks, not three — columns A–G hold *two*
vertically stacked study tables, the second titled at `A57` with its header at row 58 — and
**rows 57 and 58 are simultaneously those metadata rows and ordinary gene rows**
(`J57=KCTD5`, `Q57=PIGS`, `J58=MLST8`, `Q58=PIGU`). A header-inference reader drops those four
genes. All four are present in the TSV.

### 2.2 Two findings against the converter

1. **`wu2025_dp_overlap` merges two different variables.** Column M is *"Overlap status with
   recurrent, DP-specific, **deleted** genes"*; column T is *"…**amplified** genes"*. The
   converter writes both into one column. Nothing is lost (`role` recovers the direction) but
   the name hides it — the same shape as the `.peak_altered()` bug recorded in `CLAUDE.md`,
   where a direction column turned out to be a group label under the wrong name.
   **Proposed:** rename to `wu2025_dp_overlap_matched`, or split into
   `wu2025_dp_deleted_overlap` / `wu2025_dp_amplified_overlap`. Low risk, one file.
2. **The `stop()` regression anchor guards 519/877 but not 108/90.** A swap of columns M and T
   would leave 519/877 intact and silently invert every directional claim. The 108/90 counts
   are already computed and printed two lines above; adding them to the `stop()` condition
   closes it. **One-line change, and it is the cheap one.**
3. The header's claim that "the paper allows this" about 105 both-list genes is an
   **inference** — Wu et al. never address both-list genes at all. Reword.

Also worth recording, because it is *better* news than the converter assumes: **no gene is
`YES` in both overlap columns.** The 108 and the 90 are fully disjoint, so Wu's own directional
claim is safe on their data. The 105 both-list genes only threaten a claim built from the
*unfiltered* catalogue.

---

## 3. Why a straight transplant is wrong for ATTEND

### 3.0 ATTEND CANNOT REPRODUCE "DP-SPECIFIC" AT ALL — the design is cross-sectional

**This is the load-bearing constraint and it precedes every other consideration below.**

Every number in Wu et al.'s Figure 1 is *DP-specific*: a CNV present in the progression
tumour and **absent from that same patient's own pre-treatment baseline**. Their unit of
analysis is a **within-patient subtraction** over patient-matched trios (normal / baseline /
progression). ATTEND has no such axis. Verified three ways:

| Check | Result |
|---|---|
| Any timepoint / baseline / longitudinal concept in `code/*.R` | **None.** The only "paired" in the repo is *paired tumour-normal* for the TMB germline filter — a different meaning |
| `collapse_pid()` (`attend_harmonise.R:128`) | `group_by(.pid) \|> summarise(across(payload, num_or_first))` — **one row per patient, numeric columns averaged** |
| `add_response_class()` | a **landmark on `attend_cols$surv_time`**, not a genomic comparison |

The second is decisive: even where a patient had more than one tumour sequenced, the
integration point averages their values into a single row. The master is keyed on `pid` by
design, so any within-patient structure is gone before any report sees it. **There is no
baseline to subtract.**

Two things in the repo look adjacent and are not:

- **`response_class`** (responder / non-responder at 6 months) is a **between-patient**
  contrast on clinical outcome. Wu et al.'s is a **within-patient** contrast on genomes. They
  can say a deletion *appeared* under ICI; the strongest ATTEND equivalent is that
  non-responders *tend to carry* it. Different claim, far more confoundable.
- **Report 11's metastatic vs primary** is ATTEND vs TCGA — two cohorts, not paired patients.

**Consequence — the paper splits into two claims and only one is portable:**

| Wu et al. claim | Carried by | Portable to ATTEND? |
|---|---|---|
| Tumours **acquire** apoptotic-evasion CNVs *under* ICI pressure | the patient-matched subtraction | **No.** Requires paired samples ATTEND does not have and cannot recover from the master |
| Recurrently altered genes are **disproportionately immune-evasion genes**, direction-matched | a directional recurrent gene list + the catalogue | **Yes** — this is what §4 builds |

So the question this port can answer is **not** "what changed under treatment" but "are the
recurrent CNAs in MMRd-AS-high enriched for immune-evasion genes, in the expected direction?"
That is a legitimate extension of report 10's Q1, and it is **a different question from Wu et
al.'s**. The page must ask it in its own words rather than borrowing their framing, and must
not use the phrase "DP-specific" for anything ATTEND computes.

If a paired pre-/post-treatment axis is ever wanted, it is a **data-collection and
`build_master.R` question**, not a reporting one: it would require sampling ATTEND patients
longitudinally and giving the master a sample-level grain it deliberately does not have.

### 3.1 ≥3 of 17 does not transfer

ATTEND's MMRd-AS-high stratum is n ≈ 9. A ≥3-patient count is a 33% frequency there versus
17.6% in Wu's cohort, and neither number is a statistic. Report 10 already rejected frequency
cuts for exactly this reason — `q1_recurrence_table()` defines recurrent as **q on the
MMRd-high GISTIC run**, not a frequency cut, because aneuploidy-high tumours carry more of
everything and a bare count reads burden as recurrence.

**Do not import the ≥3 rule.** ATTEND keeps its GISTIC peak definition.

### 3.2 The two "recurrences" are not the same construct, and the page must say so

Wu's ≥3/17 cut returns **5,035 deleted genes** — roughly a quarter of the exome — with no
significance control. GISTIC's q is background-corrected. A figure that placed ATTEND's peak
genes beside Wu's list without stating the asymmetry would be making a claim neither source
supports. This belongs in the derivation note of every chunk that uses the catalogue, and the
justification in `00-methods.Rmd`.

### 3.3 An unscreened gene is not a negative — this is the repo's own rule

**This is the most important design constraint.** `CLAUDE.md` records three separate incidents
of the same failure (`pathogenic_by_patient()`, `mutation_status_long()`,
`pole_ultramutated()`): a default substituted for an absent measurement, read on the page as a
real negative result.

A gene absent from the CRISPR catalogue is **not "not a resister"** — it is **unscreened**.
The 23 pooled screens do not cover the genome uniformly. So:

- The overlap denominator must be the **screened** set, not all peak genes.
- A peak gene with no catalogue entry is **NA**, never `FALSE`.
- A companion counter (`crispr_counts()`) prints the partition — peak genes / in universe /
  screened / both-list dropped / unscreened — beside every figure, so the unscreened count is
  never implicit.

This is the same discipline as `mutation_counts()`, `promise_counts()` and
`response_counts()`, and it is non-negotiable.

### 3.4 Both-list genes must be resolved explicitly

105 symbols carry ≥2 concordant measurements as *both* resister and sensitizer. Wu et al. give
no handling. ATTEND must choose and say so.

**Proposed:** drop them from any directional test, and print the count — mirroring what
`gistic_feature_direction()` (`attend_classes.R:1549`) **already does** for genes appearing in
both amp and del peak lists ("left undirected", falls back to `|value|`). Reusing that
convention rather than inventing a second one is the `status_genes()` lesson: read the
configured source, do not restate it.

### 3.5 The intersection needs a test ATTEND supplies, and it is Family B

Wu's Venn carries no statistic. ATTEND cannot report a bare overlap count — an overlap of 12
genes means nothing without knowing what overlap chance produces.

Under `attend_scna`'s two-family rule this is **Family B (discovery), permanently**, corrected
with **BH**. The CRISPR catalogue is fixed *a priori* and externally, which is Family A shape —
but the peaks it intersects are chosen from ATTEND's own group labels, so no result can promote
it. Same logic as `perm_test_selected_panel()`. `family_adjust()` must refuse to pool it with
the 12-locus panel.

**The universe is the load-bearing choice** and must be declared in config, not at a call site:

| Candidate universe | Argument |
|---|---|
| Genes GISTIC assayed in that run (rows of `all_data_by_genes.txt`) | Matches what could have been a peak gene |
| …**intersected with** genes appearing anywhere in the catalogue (screened set) | Matches what could have been an overlap — **recommended**, and the only one consistent with §3.3 |
| All protein-coding genes | Wrong; inflates significance by counting unscreenable genes |

---

## 4. Proposed design

Nothing below is implemented.

### 4.1 Config — `attend_classes.R`

```r
attend_crispr <- list(
  file             = "crispr_t_cell_screens_wu2025.tsv",
  min_measurements = 2L,     # Wu et al. STAR Methods; reproduces 519/877 exactly
  # Direction matching is the whole method: loss of a resister and gain of a
  # sensitizer both promote immune escape. Never fold the two together.
  direction_of     = c(resister = "del", sensitizer = "amp"),
  drop_both_list   = TRUE,   # 105 genes at >=2; the paper does not address them
  universe         = "screened_assayed"   # see 3.5
)
```

### 4.2 Loader

`load_crispr_screens(min_measurements = attend_crispr$min_measurements)` — applies the cutoff
as a **declared argument**, not baked into the data file; returns `gene`, `role`,
`n_measurements`, `both_list` (logical). Returns an empty frame when the TSV is absent, so
reports stay knit-safe (the `have_cols()` philosophy).

### 4.3 Overlap + test

```r
crispr_peak_overlap(run_dir, crispr = load_crispr_screens())
```

- gene → `"amp"`/`"del"` from **that run's own** peak lists via
  `gistic_feature_direction(cnv = <run>)` — already implemented, already drops
  ambiguous-direction genes.
- intersect del-peak genes ∩ resisters, amp-peak genes ∩ sensitizers.
- universe per `attend_crispr$universe`.
- returns the overlap sets, the counts, **and the partition** for `crispr_counts()`.

`crispr_overlap_test()` — hypergeometric against the declared universe, or a permutation over
peak-gene labels if the hypergeometric's independence assumption is judged untenable for
neighbouring genes inside one peak (**open question, §6**). Family B, BH.

### 4.4 Where it renders — **decision required, §6**

Report 10 is at **exactly 800 lines = `MAX_LINES`**. A Part 4 pushes it over, exactly as the
two-cut cross-tab pushed `05-response.Rmd` over. Options in §6.

### 4.5 Tests

| Test | Pins |
|---|---|
| `test_crispr_catalogue.R` | 519/877/108/90/105 reproduce; the four boundary genes KCTD5/MLST8/PIGS/PIGU present; `min_measurements` is a real argument, not a relabel |
| `test_crispr_overlap.R` | an unscreened gene is **NA, never FALSE**; both-list genes dropped from directional output; direction matching is not symmetric (resister↔del, sensitizer↔amp) |
| `test_crispr_family.R` | `family_adjust()` refuses to pool the CRISPR overlap with the Family A panel; the universe is read from config, not a call site |

The first is a **calibration** test in the sense `test_nested_selection_permutation.R` is one:
the failure it guards (a column swap inverting every directional claim) leaves every count
intact and is invisible in a diff.

---

## 5. Honest framing on the page

This is **not** a reproduction of Wu et al.'s analysis, and the gap is larger than a
parameter choice. It is their *functional filter* applied to a different recurrence
definition, on a different tumour type, in a **cross-sectional rather than longitudinal
design**, with a test they did not compute. `00-methods.Rmd` must say:

- **ATTEND has no pre-/post-treatment axis** (§3.0). Wu et al.'s result is about what tumours
  *acquire under ICI pressure*; ATTEND's would be about what recurs *within a stratum*. The
  word "DP-specific" must not appear anywhere describing an ATTEND quantity.
- ATTEND uses GISTIC q; Wu et al. used a ≥3/17 patient count. Not the same construct.
- The overlap test is ATTEND's addition; the source reports bare Venn counts.
- The catalogue is melanoma/ICI-derived and pooled across mouse and human screens; its
  applicability to endometrial carcinoma is an assumption, not a result.
- The direction rule (loss↔resister, gain↔sensitizer) is Wu et al.'s and is carried over intact.

---

## 6. Open decisions — need your call

1. **Where does this render?** (a) new report, renumbering 11→12 and 12→13 — clean, but
   renumbering is what silently broke cross-references before (now guarded by
   `test_rmd_crossrefs.R`); (b) report 10 Part 4 with a named `OVER_CAP` entry; (c) fold into
   report 12, which already carries the immune/PD-L1 material and is thematically closest.
2. **Hypergeometric or permutation** for the overlap test? Genes inside one GISTIC peak are not
   independent, which is the hypergeometric's assumption. A permutation over peak-gene labels
   respects peak structure but is slower and needs a declared null.
3. **Which strata?** MMRd-high alone (report 10's Q1 target), or all four `scna_group` levels?
4. **Do the two converter fixes now** (§2.2: the 108/90 `stop()` anchor and the column rename)
   as a separate small commit, independent of this design?
5. **Fetch the missing supplement?** Tables S1–S3 and Figures S1–S7 were never downloaded.
   Table S1 (per-patient recurrent CNV gene lists) would let us check our reading of their
   method against their actual output — currently impossible.
6. **Given §3.0, is this worth building at all?** Honest framing of the tradeoff: without a
   pre/post axis the port answers a narrower question than the paper does, and its result
   cannot speak to treatment-driven evolution. Arguments for doing it anyway: the catalogue is
   already converted and verified, `gistic_feature_direction()` already supplies the
   direction map, and "are MMRd-AS-high peaks enriched for immune-evasion genes?" is a
   question report 10 raises and does not answer. Argument against: it adds a report (or an
   over-cap exemption) for an enrichment test on n ≈ 9, where the peak set itself turns over
   under leave-one-out — the same fragility `q1_recurrence_table()`'s `loo_retained_frac`
   column exists to expose. **Recommendation: check the LOO peak stability in report 10
   first.** If the MMRd-high peak set is not stable at n ≈ 9, a functional enrichment built on
   top of it inherits that instability and should not be built yet.

---

## 7. Provenance and limitations of the source material

The paper skill at `~/.claude/skills/wu-2025-immunity-paper` is labelled **draft**. Reason,
stated plainly: the independent second-opinion parser (pypdf 6.16.2) **cannot read Elsevier's
PDF at all** — every page raises `PdfReadError: More than one /FontFile found` because the font
descriptor for `HQHOJV+HelveticaNeue-Heavy` declares both `/FontFile` (Type 1) and `/FontFile3`
(CFF), which is invalid per spec. Reproduced directly outside the tool on 24 of 24 pages;
`strict=False` does not help. Strict verification therefore exits 2 and the package is labelled
a draft, which is the honest outcome.

Everything else verified:

- All 24 pages reviewed page-by-page against rendered previews; 0 unresolved per-page checks.
- 26 adjudications applied, 0 stale, each with a per-page reason checked against the source.
- All 8 figure crops pass the pixel check.
- Workbook: 59,260 cell coordinates compared, 0 mismatches.
- All 19 cross-page sentence joins re-audited under the builder's own semantics after a
  `reading_order` fix: 0 mis-targeted.

**Material limitations carried into every claim above:**

- **The supplement is absent.** Figures S1–S7 and Tables S1, S2, S4–S7 were never downloaded.
  The membership of the 108/90 overlap sets and the per-patient recurrent gene lists are
  therefore **unverified**.
- **In-figure text inside `figure-1`, `figure-2` and `graphical-abstract` was reviewed for panel
  completeness, not transcribed.** Venn counts quoted here were read off the delivered crop at
  magnification and are legible; finer in-panel labels are not independently verified.
- Page reviewers worked from 922×1197 previews plus targeted 2–4× crops, not the PDF at native
  resolution, so there is no character-level check of glyph variants.
