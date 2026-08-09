# data/

Read-only raw inputs. Treat everything here as immutable — never edit in place.

Expected (referenced by the reports):

- `attend_barcodes.csv` — canonical barcode → patient crosswalk (validated in report 01).

Large WES/clinical inputs live on HPC (`/hpcnfs/scratch/P_DIMA_ATTEND`) and are read by
the `code/load_*.R` loaders, not committed here. Data files are git-ignored by default
(see `.gitignore`); commit only small reference inputs.
