scripts/
--------
Analysis scripts to reproduce results. Run in numeric order (`S01` ... `S13`) after `00_check_setup.R`.

Run order
---------
- `00_check_setup.R`   Validate packages and required input files.
- `check_compadre_latest_concordance.R`  Optional concordance check against latest public COMPADRE release.
- `S01_*`              Data corrections and preprocessing (COMPADRE fixes, Ellis data).
- `S02_*`              Study selection and filtering.
- `S03_*`              Single MPM sampling distribution exploration.
- `S04_*`              Case study 1 (study-level) preprocessing and sampling caches.
- `S05_*`              Case study 1 (species-level) preprocessing and sampling caches.
- `S06_*`              Derived summaries for case study 1.
- `S07_*`              Case study 1 analysis and figures.
- `S08_*`              Case study 1 species-level supplementary analysis.
- `S09_*`              Survival issue exploration.
- `S10_*`              Variance component models.
- `S11_*`              PRISM climate download/extraction (must run before `S12_*`).
- `S12_*`              Case study 2 analysis.
- `S13_*`              Supplementary figures and summaries.
- `99_make_tables.R`   Export publication tables as DOCX.

Traceability map (manuscript claims -> code -> outputs)
--------------------------------------------------------
| Manuscript item | Script(s) | Primary outputs |
|---|---|---|
| Fig. 1 (single MPM sampling distributions) | `S03_*` | `figures/fig1_top.png`, `figures/fig1_bottom.png`; `data/derived/analysis_cache/fig1_derived_param_summary.csv` |
| Fig. 2 (shape/L and other parameter distributions) | `S07_*` | `figures/fig2_shape_l0_distributions.png`, `figures/fig2_other_param_distributions.png`; `case1_shape_summary.csv`, `case1_other_summary.csv` |
| Fig. 3 (shape-L regression and beta posterior) | `S07_*` | `figures/fig3_shape_l0_regression.png`; `data/derived/analysis_cache/case1_beta_summary.csv` |
| Table 1 (variance components) | `S10_*` | `data/derived/analysis_cache/case1_variance_ratios.csv`; `docs/tables/Table1_variance_components.docx` |
| Fig. 4 (climate effects on recruitment) | `S12_*` | `figures/clim_1.png`, `figures/clim_2.png`, `figures/clim.png`; `case2_gprc_beta_summary.csv`, `case2_spring_beta_summary.csv` |
| Target-study selection counts | `S02_*` | `data/derived/studies/target_studies.csv` |
| Supplementary Table S1-2 | `S03_*`, `99_make_tables.R` | `fig1_derived_param_summary.csv`; `docs/tables/Table_S1-2_single_mpm_derived.docx` |
| Supplementary Table S1-3 | `S13_*`, `99_make_tables.R` | `fig1_mean_mpm_param_summary.csv`; `docs/tables/Table_S1-3_mean_mpm_derived.docx` |
| Supplementary boundary figures | `S13_*` | `figures/boundary.png`, `figures/boundary2.png` |
| Supplementary species-level diagnostics | `S08_*` | `figures/sds_shape_spp.png`, `figures/sd_other_spp.png`, `figures/shape_spp.png`, `figures/shape_l0_scatter_spp.png`, `figures/hazard_trajectories_spp.png` |

Output locations
----------------
- Reused intermediate artefacts: `data/derived/analysis_cache/`
- Study-level derived data: `data/derived/studies/`
- Climate derived data: `data/derived/climate/` (notably `species_clim_prism.csv` from `S11_*`)
- Figures: `figures/`
- Publication-ready tables (DOCX): `docs/tables/`

Manuscript files
----------------
- Main manuscript: `docs/manuscript/manuscript_main.Rmd`
- Supplementary materials: `docs/manuscript/manuscript_supplement.Rmd`

Notes
-----
- Scripts load common packages/settings from `code/setup.R` (`setup_packages()`, `setup_rstan()`, `setup_prism()`).
- PRISM uses the `sf` + `terra` stack; rasters are stored under `data/raw/prism/`.
- Set `FAST_RUN=1` to shorten Stan runs in `S08_case_study_1_analysis_spp.R` during development.
