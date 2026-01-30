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
- Large downloads; store outputs outside the repo and reference in `scripts/S12_case_study_2_prism.R`.
