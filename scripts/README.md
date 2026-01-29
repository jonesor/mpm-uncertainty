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
- `S11_*`   Case study 2 analysis.
- `S12_*`   PRISM climate extraction and preprocessing.
- `S13_*`   Appendix figures/analysis.

Outputs
-------
- Intermediate artifacts are written to `data/derived/analysis_cache/`.
- Figures are written to `figures/` (many `ggsave` calls are commented).

Setup
-----
- Run `scripts/00_check_setup.R` to verify packages and required input files.
- Scripts load packages via `code/setup.R` using `setup_packages(...)`.
