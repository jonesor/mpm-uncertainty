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
| Fig. 1 single-MPM sampling distributions | `S03_*` | `figures/fig1_top.png`, `figures/fig1_bottom.png` |
| Fig. 2 parameter distributions | `S07_*` | `figures/fig2_shape_l0_distributions.png`, `figures/fig2_other_param_distributions.png` |
| Fig. 3 pace-shape regression/posteriors | `S07_*` | `figures/fig3_shape_l0_regression.png` |
| Table 1 variance components | `S10_*`, `99_make_tables.R` | `data/derived/analysis_cache/case1_variance_ratios.csv`; `docs/tables/Table1_variance_components.docx` |
| Fig. 4 climate-recruitment results | `S12_*` | `figures/clim_1.png`, `figures/clim_2.png`, `figures/clim.png` |

## Supplement traceability
| Supplement element | Script(s) | Output file(s) |
|---|---|---|
| Table S1-2 (single MPM derived parameters) | `S03_*`, `99_make_tables.R` | `data/derived/analysis_cache/fig1_derived_param_summary.csv`; `docs/tables/Table_S1-2_single_mpm_derived.docx` |
| Table S1-3 (mean MPM derived parameters) | `S13_*`, `99_make_tables.R` | `data/derived/analysis_cache/fig1_mean_mpm_param_summary.csv`; `docs/tables/Table_S1-3_mean_mpm_derived.docx` |
| Case study 1 model summary table | `S07_*` | `data/derived/analysis_cache/case1_beta_summary.csv` |
| Case study 2 model summary tables | `S12_*` | `data/derived/analysis_cache/case2_gprc_beta_summary.csv`, `data/derived/analysis_cache/case2_spring_beta_summary.csv` |
| Boundary-estimate supplementary figures | `S13_*` | `figures/boundary.png`, `figures/boundary2.png` |
| Species-level supplementary diagnostics | `S08_*` | `figures/sds_shape_spp.png`, `figures/sd_other_spp.png`, `figures/shape_spp.png`, `figures/shape_l0_scatter_spp.png`, `figures/hazard_trajectories_spp.png` |

## Supporting outputs (not directly cited as main figures/tables)
- Case study 1 diagnostics from `S06_*`, `S07_*`, `S10_*`:
  - `figures/case1_example_survival.png`
  - `figures/case1_pace_shape.png`
  - `figures/case1_varcomp_theta.png`
  - `figures/case1_varcomp_theta_summary.png`
- Case study 2 diagnostics from `S12_*`:
  - `figures/case2_temp_fecundity.png`
  - `figures/case2_beta_summary.png`

## Consistency notes
- Main manuscript now carries only **main-paper table output** (Table 1).
- Supplementary tables/figures are rendered from `docs/manuscript/manuscript_supplement.Rmd`.
- `scripts/S02_target_studies.R` still writes an extra local `studies_check.csv` artifact in the repo root; it is not used in manuscript rendering.
