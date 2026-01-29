Sources
=======

This file records dataset provenance and access notes for materials used in the paper.

COMPADRE
--------
- Dataset: COMPADRE Plant Matrix Database (pre-release)
- Build/compile date: 2017-11-22
- Notes: This pre-release dataset can be released with the paper. For long-term reproducibility, we should update analyses to the latest public COMPADRE release and record the exact version/DOI here.
- Files in repo: `data/raw/compadre/COMPADRE_v.X.X.X.RData`, `data/raw/compadre/COMPADRE_v.X.X.X_Corrected.RData`
- Script creating corrected file: `scripts/S01_compadre_correct.R`

Ellis et al. (2012)
-------------------
- Source: Ellis et al. 2012, Ecology (supplemental matrices and metadata)
- Files in repo: `data/raw/ellis_2012/Transition_Matrices.txt`, `Population_data.txt`, `Species_Information.txt`, `metadata.htm`, `default.htm`

PRISM climate data
------------------
- Source: PRISM Climate Group, Oregon State University (monthly tmean and ppt rasters)
- Download helpers: `scripts/download/fetch_prism_tmp.sh`, `scripts/download/fetch_prism_ppt.sh`
- Derived outputs: `data/derived/climate/species_clim_prism.csv`

Study-level inputs
------------------
- Source list and DOIs: `data/derived/studies/_data_sources.csv`
- Derived inputs: `data/derived/studies/*_n.csv`, `data/derived/studies/*_sim.RData`
