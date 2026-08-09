# output/

Processed, regenerable products of the pipeline.

- `clean_data/attend_master_joined.csv` — the per-patient master table written by
  report 02 and read by reports 03–06. This is the integration boundary of the pipeline.
- `clean_data/<dataset>.csv` — each raw table written out individually by report 02.

Contents are git-ignored by default because they are reproducible from the reports.
