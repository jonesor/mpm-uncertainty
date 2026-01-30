
### libraries
source("code/setup.R")
setup_packages(c("tidyverse", "sf", "terra", "Rcompadre", "prism"))
source("code/functions.R")


### create shell scripts to download prism climate rasters
# year <- 1981:2013
# base_ppt <- "wget ftp://prism.nacse.org/monthly/ppt/"
# base_tmp <- "wget ftp://prism.nacse.org/monthly/tmean/"
# mid_ppt <- "/PRISM_ppt_stable_4kmM3_"
# mid_tmp <- "/PRISM_tmean_stable_4kmM2_"
# end <- "_all_bil.zip"
# 
# cat(c("#!/usr/bin/env bash", paste0(base_ppt, year, mid_ppt, year, end)),
#     file = "scripts/download/fetch_prism_ppt.sh", sep = "\n")
# cat(c("#!/usr/bin/env bash", paste0(base_tmp, year, mid_tmp, year, end)),
#     file = "scripts/download/fetch_prism_tmp.sh", sep = "\n")


### prism raster files
if (!dir.exists("prism")) {
  dir.create("prism", recursive = TRUE)
}

prism::prism_set_dl_dir("prism")

maybe_download_prism <- function(years) {
  if (length(list.files("prism", pattern = "\\.bil$")) > 0) {
    return(invisible(NULL))
  }
  prism::get_prism_monthlys(type = "ppt", years = years, mon = 1:12, keepZip = FALSE)
  prism::get_prism_monthlys(type = "tmean", years = years, mon = 1:12, keepZip = FALSE)
}

prism_files <- paste0("prism/", list.files("prism"))
month_files <- prism_files[grepl("[[:digit:]]{6}", prism_files)]
bil_files <- month_files[grepl(".bil$", month_files)]
files_ppt <- bil_files[grepl("ppt", bil_files)]
files_tmp <- bil_files[grepl("tmean", bil_files)]


### coordinates for each species/pop of interest
coords_path <- "data/derived/climate/species_coords.csv"
if (!file.exists(coords_path)) {
  compadre <- cdb_fetch("data/raw/compadre/COMPADRE_v.X.X.X_Corrected.RData")
  comp_sub <- compadre %>% 
    filter(MatrixComposite == "Individual",
           MatrixTreatment == "Unmanipulated",
           ProjectionInterval == "1",
           MatrixCaptivity == "W")
  comp_time_series <- comp_sub %>% 
    as_tibble() %>% 
    filter(!is.na(Lon) & !is.na(Lat)) %>% 
    group_by(SpeciesAuthor) %>% 
    mutate(n_year = length(unique(MatrixStartYear))) %>% 
    ungroup() %>% 
    filter(n_year >= 5) %>% 
    group_by(SpeciesAuthor, MatrixPopulation) %>%
    summarize(Lon = unique(Lon)[1],
              Lat = unique(Lat)[1],
              n_year = unique(n_year),
              .groups = "drop")
  write_csv(comp_time_series, coords_path)
}

spp_df <- read_csv(coords_path) %>% 
  filter(SpeciesAuthor == "Silene_spaldingii")

if (length(prism_files) == 0) {
  compadre <- cdb_fetch("data/raw/compadre/COMPADRE_v.X.X.X_Corrected.RData")
  years <- compadre %>% 
    as_tibble() %>% 
    filter(SpeciesAuthor == "Silene_spaldingii") %>% 
    filter(!is.na(MatrixStartYear)) %>% 
    pull(MatrixStartYear) %>% 
    unique() %>% 
    sort()
  maybe_download_prism(years)
  prism_files <- paste0("prism/", list.files("prism"))
}


### get climate data from all raster files for all species of interest
df_clim <- tibble(file_tmp = files_tmp, file_ppt = files_ppt) %>% 
  mutate(Date = map_chr(files_tmp, ~ strsplit(.x, "_")[[1]][5])) %>% 
  mutate(Year = as.integer(substr(Date, 1, 4))) %>% 
  mutate(Month = as.integer(substr(Date, 5,6))) %>% 
  group_by(Year, Month) %>% 
  do(fetch_prism(.$file_tmp, .$file_ppt, spp_df)) %>% 
  ungroup() %>% 
  arrange(SpeciesAuthor, MatrixPopulation, Year, Month) %>% 
  filter(!(is.na(tmp) & is.na(ppt)))

write.csv(df_clim, "data/derived/climate/species_clim_prism.csv", row.names = FALSE)
