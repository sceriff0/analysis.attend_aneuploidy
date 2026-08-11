# Runbook — classification options beyond the current GISTIC reproduction

**Status:** backlog, not implemented. **Relates to:** `analysis/40_tcga_classification.Rmd`.

Where the copy-number classification could go beyond reproducing Kandoth et al. Suppl.
Methods S2, given the tools this project already runs (DRAGEN, ASCETS, variantalker).
Kept out of the report because none of it changes what that report currently shows.


Ordered by effort-to-payoff:

1. **GISTIC thresholded-peak clustering — DONE (the only clustering).** Reproduces Kandoth
   Suppl. S2 exactly. The **ASCETS** aneuploidy score (`data/ascets/`) is still consumed as
   an **independent** per-sample aneuploidy value for the enrichment cross-check; it does
   not enter the clustering.
2. **Metric/linkage — this variant uses 1 − Pearson / Ward.D2.** S2 leaves the CN
   metric/linkage unstated; report 40 applies the mRNA/methylation pair (`cnv_dist =
   "correlation"`, `cnv_method = "ward.D2"`). The k = 2..8 sweep and the
   silhouette-vs-k panel make the cut choice inspectable; compare against base 09 (Euclidean/complete) for robustness.
3. **POLE branch removed.** ATTEND has no POLE-ultramutated cases, so the cascade now
   starts at MSI. `pole_ultramutated()` / `attend_pole` are kept (unwired) — if a POLE
   branch is ever wanted, feed a MAF with a protein-change column to `pole_ultramutated()`
   (matches P286R/V411L directly) and reinstate the branch in `add_tcga_class()`.
4. **Purity/ploidy-correct the copy number (medium effort).** If DRAGEN ran
   tumor-only, log2 ratios aren't purity-corrected, blurring low-purity samples.
   Re-call with **ASCAT/FACETS** (matched normal) for integer allele-specific CN,
   whole-genome-doubling and LOH — sharpening the CN-high boundary and enabling a
   real **HRD** cross-check against `hrd__HRD_Score`.
5. **Serous-like = SCNA cluster membership (faithful TCGA), TP53 validated not gated.**
   `add_tcga_class()` sets the CN-high call via `attend_tcga_serous_rule`, now defaulting
   to **`"cluster_only"`** — serous-like = the high-aneuploidy CNV cluster, exactly as TCGA
   defined it (Fig. 2b, cluster 4). *TP53* is checked as a **property** (the TP53-enrichment
   QC in report 40 should show ~90% mutant), not required. Alternatives remain: `"and"`
   (also require TP53, higher precision), `"or"`, `"tp53_only"`. The top of the cascade is
   **MMRd = MMR-loss by IHC (MMRd vs MMRp)**; MSI-high (NGS) is a correlate, not a criterion.
   Remaining focal-CNV upgrade is note 2 (GISTIC).
6. **Ancestry-aware TMB — robust/matched DONE; `nassar` scaffolded.** The recompute
   (`tmb_variant_table()`) produces ancestry-robust (grpmax) / ancestry-matched (on the
   upstream TMB_SCORE baseline)
   TMB from the gnomAD-annotated MAFs (matched activates once
   `attend_ancestry$clinical_col` points at the real ancestry field; swap in
   EthSEQ/somalier genetic ancestry later for accuracy). The **`nassar`** regression
   variant (`apply_nassar_recalibration()`) is wired with identity coefficients — the
   remaining step is to **fit `(m, b)` per ancestry on TCGA WES** (paired vs
   tumor-only-simulated TMB) and paste them into `attend_tmb_nassar$coef`.

