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
- `S11_*`   PRISM climate extraction and preprocessing (includes `prism` package download workflow; run before `S12_*`).
- `S12_*`   Case study 2 analysis (run after `S11_*` produces climate inputs).
- `S13_*`   Appendix figures/analysis.
- `99_*`    Manuscript table exports (DOCX) via `officer` + `flextable`.

Outputs
-------
- Intermediate artifacts are written to `data/derived/analysis_cache/`.
- Figures are written to `figures/` with filenames set in the scripts.
- Climate extraction outputs are written to `data/derived/climate/` (notably `species_clim_prism.csv` from `S11_*`).

Figures and tables
------------------
- **Fig. 1 (sampling distributions for a single MPM)**: `S03_*` (`figures/fig1_top.png`, `figures/fig1_bottom.png`).
- **Fig. 2 (sampling distributions vs point estimates for shape/life expectancy)**: `S07_*` (`figures/fig2_shape_l0_distributions.png`).
- **Fig. 2 (other parameter distributions)**: `S07_*` (`figures/fig2_other_param_distributions.png`).
- **Fig. 3 (shape–life expectancy relationship + beta posterior)**: `S07_*` (`figures/fig3_shape_l0_regression.png`).
- **Table 1 (variance components)**: `S10_*` (exported to `data/derived/analysis_cache/case1_variance_ratios.csv`).
- **Fig. 4 (Silene climate analyses)**: `S12_*` (`figures/clim_1.png`, `figures/clim_2.png`, `figures/clim.png`; plus `case2_*` diagnostic figures).
- **Appendix figures (boundary estimates)**: `S13_*` (`figures/boundary.png`, `figures/boundary2.png`).
- **Species‑level case study figures**: `S08_*` (`figures/sds_shape_spp.png`, `figures/sd_other_spp.png`, `figures/shape_spp.png`, `figures/shape_l0_scatter_spp.png`, `figures/hazard_trajectories_spp.png`) for supplementary material.

Setup
-----
- Run `scripts/00_check_setup.R` to verify packages and required input files.
- Scripts load packages via `code/setup.R` using `setup_packages(...)`.
- PRISM download location is configured by `setup_prism()` in `code/setup.R`.
- PRISM extraction uses the `sf` + `terra` stack (no `rgdal`).
- Set `FAST_RUN=1` in the environment to shorten Stan runs in `S08_case_study_1_analysis_spp.R` for development.

Linting and formatting
----------------------
- Linting uses `.lintr` in the repo root (object-usage checks are disabled; naming and pipe preferences are not enforced).
- Run from R: `lintr::lint_dir(".")`.
- Formatting uses `scripts/format.R` (wrapper around `styler`).
