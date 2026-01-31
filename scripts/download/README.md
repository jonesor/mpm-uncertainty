scripts/download/
-----------------
Helper scripts to download external climate data (PRISM). These are not run automatically.

Contents
--------
- `fetch_prism_ppt.sh`  Download monthly precipitation rasters.
- `fetch_prism_tmp.sh`  Download monthly mean temperature rasters.

Notes
-----
- Scripts use FTP URLs; access may be blocked. Prefer the `prism` R package workflow embedded in `scripts/S12_case_study_2_prism.R`.
- Downloads are large; by default `prism` writes to `data/raw/prism/`, but you can set another download directory in `scripts/S12_case_study_2_prism.R`.
