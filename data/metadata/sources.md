Sources
=======

This file records dataset provenance and access notes for materials used in the paper.

COMPADRE
--------
- Dataset: COMPADRE Plant Matrix Database
- Version used in analyses: `6.26.3.0`
- Release date: `Mar_11_2026`
- Files in repo: `data/raw/compadre/COMPADRE_v6.26.3.0.RData`, `data/raw/compadre/COMPADRE_v6.26.3.0_Corrected.RData`
- Script creating corrected file: `scripts/S01_compadre_correct.R`
- Notes: The corrected file applies project-level fixes documented in `scripts/S01_compadre_correct.R`.

Ellis et al. (2012)
-------------------
- Source: Ellis et al. 2012, Ecology (supplemental matrices and metadata)
- Files in repo: `data/raw/ellis_2012/Transition_Matrices.txt`, `Population_data.txt`, `Species_Information.txt`

PRISM climate data
------------------
- Source: PRISM Climate Group, Oregon State University (monthly tmean and ppt rasters)
- Download helper: `scripts/S11_case_study_2_prism.R` via `setup_prism()` in `code/setup.R` for analysis 2
- Derived outputs: `data/derived/climate/species_clim_prism.csv`
- Processing stack: uses `sf` + `terra`.

Study-level inputs
------------------
- Source list and DOIs: `data/derived/studies/_data_sources.csv`
- Derived inputs: `data/derived/studies/*_n.csv`, `data/derived/studies/*_sim.RData`
