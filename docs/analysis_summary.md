# Analysis Traceability Summary

This file maps manuscript statements and display items to the scripts and concrete output files that generate them.

## Scope and manuscript split
- **Main manuscript:** `docs/manuscript/manuscript_main.Rmd`
- **Supplementary materials:** `docs/manuscript/manuscript_supplement.Rmd`
- **Pipeline scripts:** `scripts/S01_*` ... `scripts/S13_*` (+ `scripts/00_check_setup.R`, `scripts/99_make_tables.R`)

## End-to-end order
1. `scripts/00_check_setup.R`
2. `scripts/S01_*` ... `scripts/S10_*`
3. `scripts/S11_*` (PRISM extraction)
4. `scripts/S12_*` (climate analysis)
5. `scripts/S13_*` (supplement summaries)
6. `scripts/99_make_tables.R` (DOCX tables)

## Main manuscript traceability
| Manuscript element | Script(s) | Output file(s) |
|---|---|---|
| Data corrections and COMPADRE preprocessing (Methods) | `S01_*` | `data/raw/compadre/COMPADRE_v.X.X.X_Corrected.RData` |
| Target-study filtering (Methods) | `S02_*` | `data/derived/studies/target_studies.csv` |
| Fig. 1 single-MPM sampling distributions | `S03_*` | `figures/Figure_1_transition_rates_single_mpm.png`, `figures/Figure_2_derived_parameters_single_mpm.png` |
| Fig. 2 parameter distributions | `S07_*` | `figures/Figure_3_analysis1_point_vs_sampling_distributions.png`, `figures/Analysis1_additional_parameter_distributions.png` |
| Fig. 3 pace-shape regression/posteriors | `S07_*` | `figures/Figure_4_analysis1_life_expectancy_shape_relationship.png` |
| Table 1 variance components | `S10_*`, `99_make_tables.R` | `data/derived/analysis_cache/case1_variance_ratios.csv`; `docs/tables/Table1_variance_components.docx` |
| Fig. 4 climate-recruitment results | `S12_*` | `figures/Analysis2_recruitment_model_fits.png`, `figures/Analysis2_monthly_lag_coefficients.png`, `figures/Figure_5_analysis2_climate_effects_recruitment.png` |

## Supplement traceability
| Supplement element | Script(s) | Output file(s) |
|---|---|---|
| Table S1-2 (single MPM derived parameters) | `S03_*`, `99_make_tables.R` | `data/derived/analysis_cache/fig1_derived_param_summary.csv`; `docs/tables/Table_S1-2_single_mpm_derived.docx` |
| Table S1-3 (mean MPM derived parameters) | `S13_*`, `99_make_tables.R` | `data/derived/analysis_cache/fig1_mean_mpm_param_summary.csv`; `docs/tables/Table_S1-3_mean_mpm_derived.docx` |
| Case study 1 model summary table | `S07_*` | `data/derived/analysis_cache/case1_beta_summary.csv` |
| Case study 2 model summary tables | `S12_*` | `data/derived/analysis_cache/case2_gprc_beta_summary.csv`, `data/derived/analysis_cache/case2_spring_beta_summary.csv` |
| Boundary-estimate supplementary figures | `S13_*` | `figures/Figure_S1_boundary_estimate_diagnostic.png`, `figures/Figure_S2_boundary_survivorship_illustration.png` |
| Species-level supplementary diagnostics | `S08_*` | `figures/Analysis1_species_shape_life_expectancy_distributions.png`, `figures/Analysis1_species_other_parameter_distributions.png`, `figures/Analysis1_species_shape_life_expectancy_regression.png`, `figures/Analysis1_species_shape_life_expectancy_scatter.png`, `figures/Analysis1_species_hazard_trajectories.png` |

## Supporting outputs (not directly cited as main figures/tables)
- Case study 1 diagnostics from `S06_*`, `S07_*`, `S10_*`:
  - `figures/Analysis1_example_survivorship_and_hazard.png`
  - `figures/Analysis1_pace_shape_scatter.png`
  - `figures/Analysis1_variance_components_scatter.png`
  - `figures/Analysis1_variance_components_summary.png`
- Case study 2 diagnostics from `S12_*`:
  - `figures/Analysis2_spring_temperature_recruitment_scatter.png`
  - `figures/Analysis2_spring_temperature_beta_summary.png`

## Consistency notes
- Main manuscript now carries only **main-paper table output** (Table 1).
- Supplementary tables/figures are rendered from `docs/manuscript/manuscript_supplement.Rmd`.
- `scripts/S02_target_studies.R` still writes an extra local `studies_check.csv` artifact in the repo root; it is not used in manuscript rendering.
