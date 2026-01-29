### preflight checks for dependencies and required files
source("code/setup.R")

required_pkgs <- c(
  "tidyverse",
  "Rcompadre",
  "Rage",
  "popbio",
  "popdemo",
  "gridExtra",
  "cowplot",
  "ggridges",
  "rstan",
  "loo",
  "sf",
  "terra"
)

setup_packages(required_pkgs)
setup_rstan()

required_files <- c(
  "data/raw/compadre/COMPADRE_v.X.X.X.RData",
  "data/raw/compadre/COMPADRE_v.X.X.X_Corrected.RData",
  "data/raw/ellis_2012/Transition_Matrices.txt",
  "data/derived/studies/_data_sources.csv",
  "data/derived/studies/aschero_U.RData",
  "data/derived/climate/species_coords.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste("-", missing_files, collapse = "\n"),
    call. = FALSE
  )
}

message("Preflight checks passed.")
