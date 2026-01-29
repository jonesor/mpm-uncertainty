
### libraries
source("code/setup.R")
setup_packages(c("tidyverse", "rgdal", "raster", "sp"))


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
prism_files <- paste0("prism/", list.files("prism"))
month_files <- prism_files[grepl("[[:digit:]]{6}", prism_files)]
bil_files <- month_files[grepl(".bil$", month_files)]
files_ppt <- bil_files[grepl("ppt", bil_files)]
files_tmp <- bil_files[grepl("tmean", bil_files)]


### coordinates for each species/pop of interest
spp_df <- read_csv("data/derived/climate/species_coords.csv") %>% 
  filter(SpeciesAuthor == "Silene_spaldingii")


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
