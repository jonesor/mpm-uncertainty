data/derived/climate/
---------------------
Derived climate datasets used in case study 2.

Contents
--------
- `species_coords.csv`       Species/population coordinates for climate extraction.
- `species_clim_prism.csv`   Extracted PRISM temperature/precipitation summaries.

Notes
-----
PRISM raster downloads and extraction are handled by `scripts/S12_case_study_2_prism.R` using the `prism` package (rasters in `data/raw/prism/`); extraction uses `sf` + `terra`.
