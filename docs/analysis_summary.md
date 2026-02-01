# Analysis Summary (manuscript scaffold)

This scaffold aligns the manuscript sections with the analysis scripts (S01–S13) and highlights where figures, tables, and model summaries are produced. Add concrete results as they are finalized; placeholders below mark expected outputs.
Numeric exports are rounded to manuscript‑appropriate precision (typically 3 significant figures; probabilities to 3 decimals).

## Title/Abstract
- **Scripts:** S01–S13 (overall pipeline).
- **Add later:** One‑paragraph summary of main findings; include key effect sizes and uncertainty statements.

## 1. Introduction
- **Primary sources:** manuscript text only.
- **No script outputs** (contextual framing).

## 2. Methods

### 2.1 COMPADRE data and corrections
- **Scripts:** `S01_*` (COMPADRE fixes), `S02_*` (study selection).
- **Possible outputs to export:** summary of corrections (CSV/MD), final selection counts by filter step.

### 2.2 Target studies
- **Scripts:** `S02_*`
- **Exports:** `data/derived/studies/target_studies.csv`.

### 2.3 Typical construction of MPMs
- **Scripts:** narrative only; data example in `S03_*`.
- **Possible outputs:** example transition counts used in Fig. 1 (CSV/MD).

### 2.4 Sampling distribution for a single MPM
- **Scripts:** `S03_*`
- **Figures:** **Fig. 1** (`figures/fig1_top.png`, `figures/fig1_bottom.png`).
- **Exports:** `data/derived/analysis_cache/fig1_derived_param_summary.csv`.

### 2.5 Case Study #1: mortality trajectories
- **Scripts:** `S04_*`, `S05_*`, `S06_*`, `S07_*`, `S10_*`, `S13_*`
- **Figures:** 
  - **Fig. 2** (`figures/fig2_shape_l0_distributions.png`, `figures/fig2_other_param_distributions.png`)
  - **Fig. 3** (`figures/fig3_shape_l0_regression.png`)
- **Tables:** 
  - **Table 1** variance components (`S10_*`; exported to `data/derived/analysis_cache/case1_variance_ratios.csv`).
- **Supplement:** boundary estimates and related diagnostics (`S13_*`).

### 2.6 Case Study #2: weather impacts
- **Scripts:** `S11_*` (PRISM extraction), `S12_*` (analysis)
- **Figures:** **Fig. 4** (`S12_*`: `figures/clim_1.png`, `figures/clim_2.png`, `figures/clim.png`).
- **Exports:** `data/derived/analysis_cache/case2_gprc_beta_summary.csv`, `data/derived/analysis_cache/case2_spring_beta_summary.csv`.

## 3. Results

### 3.1 Comparative analysis of mortality trajectories
- **Scripts:** `S06_*`, `S07_*`, `S10_*`
- **Include:** key effect estimates, posterior probabilities, variance components.
- **Recommended outputs:** 
  - `case1_variance_ratios.csv` (from `S10_*`)
  - `case1_beta_summary.csv` (from `S07_*`)

### 3.2 Case Study #2: weather and recruitment
- **Scripts:** `S12_*`
- **Include:** GPRC lag effects, simplified Spring temperature model results.
- **Recommended outputs:** 
  - `case2_gprc_beta_summary.csv`
  - `case2_spring_beta_summary.csv`

## 4. Discussion
- **Scripts:** narrative only.
- **Optionally cite:** robustness checks in `S13_*` and sensitivity implied by variance components.

## Data Accessibility
- **Scripts:** `S01_*` (corrected COMPADRE), `S11_*` (PRISM extraction).
- **Add later:** repository + archive DOI once finalized.

## Appendices

### Appendix S1
- **Scripts:** `S03_*` (single MPM sampling), `S13_*` (boundary survival).
- **Figures:** `figures/boundary.png`, `figures/boundary2.png`.

### Appendix S2
- **Scripts:** `S07_*`, `S08_*`, `S10_*`, `S12_*` (Stan model details and diagnostics).
- **Recommended outputs:** 
  - `table_mcmc_settings.md` (model settings and priors)
  - `table_model_spec_case1.md`, `table_model_spec_case2.md`

## Suggested additional exports
- **S02:** CSV table of included studies and sample sizes (now `data/derived/studies/target_studies.csv`).
- **S06/S07:** CSV of point estimates vs sampling distribution medians/CI (now `case1_shape_summary.csv`, `case1_other_summary.csv`, `case1_beta_summary.csv`).
- **S10:** CSV of variance components (now `case1_variance_ratios.csv`).
- **S12:** CSV of lag coefficients and posterior summaries for Fig. 4 (now `case2_gprc_beta_summary.csv`, `case2_spring_beta_summary.csv`).
