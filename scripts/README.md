scripts/
--------
Analysis scripts to reproduce results. Intended to be run in numeric order (S01 ... S13).

Contents
--------
- `S01_*`   Data corrections and preprocessing (COMPADRE fixes, Ellis data).
- `S02_*`   Study selection and filtering.
- `S03_*`   Single MPM sampling distribution exploration.
- `S04_*`   Case study 1 (study-level) preprocessing and caching.
- `S05_*`   Case study 1 (species-level) preprocessing and caching.
- `S06_*`   Derived summaries for case study 1.
- `S07_*`   Case study 1 analysis and figures.
- `S08_*`   Case study 1 species analysis.
- `S09_*`   Survival issue exploration.
- `S10_*`   Variance component models.
- `S11_*`   Case study 2 analysis (run after `S12_*` produces climate inputs).
- `S12_*`   PRISM climate extraction and preprocessing (includes `prism` package download workflow; run before `S11_*`).
- `S13_*`   Appendix figures/analysis.

Outputs
-------
- Intermediate artifacts are written to `data/derived/analysis_cache/`.
- Figures are written to `figures/` with filenames set in the scripts.

Setup
-----
- Run `scripts/00_check_setup.R` to verify packages and required input files.
- Scripts load packages via `code/setup.R` using `setup_packages(...)`.
- PRISM extraction uses the `sf` + `terra` stack (no `rgdal`).
- Set `FAST_RUN=1` in the environment to shorten Stan runs in `S08_case_study_1_analysis_spp.R` for development.
