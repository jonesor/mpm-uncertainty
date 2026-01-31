data/derived/climate/
---------------------
Derived climate datasets used in case study 2.

Contents
--------
- `species_coords.csv`       Species/population coordinates for climate extraction.
- `species_clim_prism.csv`   Extracted PRISM temperature/precipitation summaries.

Notes
-----
PRISM raster downloads and extraction are handled by `scripts/S11_case_study_2_prism.R` using `setup_prism()` in `code/setup.R` (rasters in `data/raw/prism/`); extraction uses `sf` + `terra`.
