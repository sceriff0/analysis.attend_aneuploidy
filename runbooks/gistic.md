# Runbook — producing GISTIC output for the TCGA cascade

**Reads:** DRAGEN `.seg` files on HPC. **Writes:** `data/gistic/`.
**Consumed by:** `analysis/08_tcga_classification.Rmd`.

This is an operational recipe, run once per cohort freeze, not per knit. The rationale for
using GISTIC thresholded peaks as the clustering feature is in `analysis/00_methods.Rmd`.


The copy-number clustering in `08_tcga_classification.Rmd` **requires GISTIC output in
`data/gistic/`** and has no fallback. GISTIC is a compiled MATLAB binary that runs
**outside R**, so you produce its output once, then knit. The recipe reproduces Kandoth et
al. Suppl. Methods S2: SNP6 to tangent-normalisation to CBS to GISTIC 2.0 there, DRAGEN
`.seg` to GISTIC 2.0 here.

**Step 1 — build the GISTIC input (R, on the cluster where the `.seg` live):**

```bash
Rscript -e 'source("code/load_wes_results.R"); write_gistic_seg()'
# -> output/gistic_input/attend_all_segments.seg  (Sample/Chr/Start/End/Num_Markers/Seg.CN)
```

**Step 2 — run GISTIC via `code/run_gistic.sh`** (runs both locally and under SLURM):

```bash
# get GISTIC (pick what your HPC has):
singularity pull gistic2.sif docker://genepattern/docker-gistic   # or: module load gistic/2.0.23
# grab the hg38 refgene .mat (stock binary ships hg19 only — using it MISPLACES every peak):
#   hg38.UCSC.add_mir.160920.refgene.mat  (GISTIC2 GitHub / bzhanglab/GISTIC2_example)

# submit to the scheduler (or run ./code/run_gistic.sh interactively):
GISTIC_SIF=gistic2.sif REFGENE=/path/hg38.UCSC.add_mir.160920.refgene.mat \
  sbatch code/run_gistic.sh
```

The script carries `#SBATCH` directives (4 h, 32 GB, 4 cores), auto-wraps a Singularity
image when `GISTIC_SIF` is set (else uses a `module`/conda `gistic2`), `cd`s to the repo
root, and checks the seg + refgene exist. It writes `data/gistic/`:
`all_thresholded.by_genes.txt` (the clustering input), `amp_genes`/`del_genes.conf_*.txt`
(the significant peaks), `all_lesions.conf_*.txt` (the report-40 oncoplot overlay).

**Step 3 — knit `analysis/08_tcga_classification.Rmd`.** `load_gistic_thresholded()`
reads `all_thresholded.by_genes.txt`, restricts to the significant-peak genes, and the CNV
sections populate automatically.


