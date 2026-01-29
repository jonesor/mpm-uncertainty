MPM Sampling Uncertainty
========================

Data and code for analyses in Barks and Jones: Accounting for sampling uncertainty in analyses of published projection matrices.

This repository is organized to support replicability and FAIR principles. The main analysis workflow lives in `scripts/`, reusable functions in `code/`, models in `models/`, and data split into `data/raw/` (as obtained) and `data/derived/` (processed/intermediate).

Quick start
-----------
- Run `scripts/00_check_setup.R` to verify packages and required inputs.
- Start with `scripts/` in numeric order (S01 ... S13). Each script is intended to be run after prior steps.
- Outputs that are reused across scripts are cached under `data/derived/analysis_cache/`.
- Figures are written to `figures/` (commented `ggsave` calls indicate intended outputs).

Data sources
------------
- COMPADRE: pre-release build compiled 2017-11-22 (recorded in `data/metadata/sources.md`).
- Ellis et al. (2012) supplemental matrices.
- PRISM climate rasters (download helpers in `scripts/download/`).

Folder guide
------------
- `code/`               Reusable R functions used by multiple scripts.
- `scripts/`            Analysis scripts (numbered, run in order).
- `models/`             Stan model definitions.
- `data/raw/`           Raw source data, as obtained from external sources.
- `data/derived/`       Processed data and intermediate analysis artifacts.
- `data/metadata/`      Data dictionaries and provenance notes.
- `figures/`            Generated figures for the paper.
- `supplement/`         Non-reproducible source assets (e.g., figure source files).
- `docs/`               Notes, methods, and documentation.

Notes
-----
- Some external data sources require download (e.g., PRISM rasters). See `scripts/download/` and `data/metadata/`.
- The project uses R and Stan; package versions are not yet pinned.
