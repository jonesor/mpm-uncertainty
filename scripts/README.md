scripts/
--------
Analysis scripts to reproduce results. Run in numeric order (`S01` ... `S15`) after `00_check_setup.R`.

Run order
---------
- `00_check_setup.R`   Validate packages and required input files.
- `S01_*`              Data corrections and preprocessing (COMPADRE fixes, Ellis data).
- `S02_*`              Study selection and filtering.
- `S03_*`              Single MPM sampling distribution exploration.
- `S04_*`              Analysis 1 (study-level) preprocessing and sampling caches.
- `S05_*`              Analysis 1 (species-level) preprocessing and sampling caches.
- `S06_*`              Derived summaries for analysis 1.
- `S07_*`              Analysis 1 figures and model outputs.
- `S08_*`              Analysis 1 species-level supplementary analysis.
- `S09_*`              Survival issue exploration.
- `S10_*`              Variance component models.
- `S11_*`              PRISM climate download/extraction (must run before `S12_*`).
- `S12_*`              Analysis 2.
- `S13_*`              Supplementary figures and summaries.
- `S14_*`              Analysis 3 (Astragalus, all COMPADRE sites) climate-recruitment analysis.
- `S15_*`              Analysis 1 sensitivity to the quasi-stable truncation threshold for shape.
- `99_make_tables.R`   Export publication tables as DOCX.

Traceability map (manuscript claims -> code -> outputs)
--------------------------------------------------------
| Manuscript item | Script(s) | Primary outputs |
|---|---|---|
| Fig. 1 (single MPM sampling distributions) | `S03_*` | `figures/Figure_1_transition_rates_single_mpm.png`, `figures/Figure_2_derived_parameters_single_mpm.png`; `data/derived/analysis_cache/fig1_derived_param_summary.csv` |
| Fig. 3 (shape/L and other parameter distributions) | `S07_*` | `figures/Figure_3_analysis1_point_vs_sampling_distributions.png`, `figures/Analysis1_additional_parameter_distributions.png`; `case1_shape_summary.csv`, `case1_other_summary.csv` |
| Fig. 4 (shape-L regression and beta posterior) | `S07_*` | `figures/Figure_4_analysis1_life_expectancy_shape_relationship.png`; `data/derived/analysis_cache/case1_beta_summary.csv` |
| Table 1 (variance components) | `S10_*` | `data/derived/analysis_cache/case1_variance_ratios.csv`; `docs/tables/Table1_variance_components.docx` |
| Fig. 5 (climate effects on recruitment + stage-specific survival) | `S12_*` | `figures/Analysis2_recruitment_model_fits.png`, `figures/Analysis2_monthly_lag_coefficients.png`, `figures/Analysis2_stage_specific_survival_beta_summary.png`, `figures/Figure_5_analysis2_climate_effects_recruitment.png`; `case2_gprc_beta_summary.csv`, `case2_spring_beta_summary.csv` |
| Fig. 6 (Astragalus multisite climate analysis) | `S14_*` | `figures/Analysis3_<site>_spring_temperature_recruitment_scatter.png`, `figures/Analysis3_<site>_spring_temperature_beta_summary.png`, `figures/Analysis3_<site>_spring_temperature_model_fits.png`, `figures/Analysis3_<site>_survival_model_fits.png`, `figures/Analysis3_<site>_survival_stage_beta_summary.png`, `figures/Analysis3_<site>_survival_stage_curves.png`, `figures/Analysis3_<site>_monthly_lag_coefficients.png`, `figures/Figure_6_analysis3_climate_effects_multisite.png`; `case3_astragalus_coverage_<site>.csv`, `case3_astragalus_spring_beta_summary.csv`, `case3_astragalus_gprc_beta_summary.csv`, `case3_astragalus_stan_diagnostics.csv`, `case3_astragalus_site_summary.csv` |
| Target-study selection counts | `S02_*` | `data/derived/studies/target_studies.csv` |
| Supplementary Table S1-2 | `S03_*`, `99_make_tables.R` | `fig1_derived_param_summary.csv`; `docs/tables/Table_S1-2_single_mpm_derived.docx` |
| Supplementary Table S1-3 | `S13_*`, `99_make_tables.R` | `fig1_mean_mpm_param_summary.csv`; `docs/tables/Table_S1-3_mean_mpm_derived.docx` |
| Supplementary boundary figures | `S13_*` | `figures/Figure_S1_boundary_estimate_diagnostic.png`, `figures/Figure_S2_boundary_survivorship_illustration.png` |
| Supplementary truncation-threshold sensitivity | `S15_*`, `99_make_tables.R` | `figures/Figure_S4_shape_threshold_sensitivity.png`; `data/derived/analysis_cache/case1_qsd_threshold_sensitivity_summary.csv`; `docs/tables/Table_S14_qsd_threshold_sensitivity.docx` |
| Supplementary species-level diagnostics | `S08_*` | `figures/Analysis1_species_shape_life_expectancy_distributions.png`, `figures/Analysis1_species_other_parameter_distributions.png`, `figures/Analysis1_species_shape_life_expectancy_regression.png`, `figures/Analysis1_species_shape_life_expectancy_scatter.png`, `figures/Analysis1_species_hazard_trajectories.png` |

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
- Plot styling is centralized in `code/functions.R` (`theme_mpm()`, `mpm_colors()`, `set_mpm_plot_defaults()`).
- PRISM uses the `sf` + `terra` stack; rasters are stored under `data/raw/prism/`.
- Set `FAST_RUN=1` to shorten Stan runs in `S08_case_study_1_analysis_spp.R` during development.
- Archived comparison scripts live in `scripts/archive/` and are not part of the canonical workflow.
