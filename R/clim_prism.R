
# libraries
library(rgdal)
library(raster)
library(readr)
library(tibble)
library(dplyr)
library(purrr)


## create shell scripts to download prism climate rasters
# year <- 1981:2013
# base_ppt <- "wget ftp://prism.nacse.org/monthly/ppt/"
# base_tmp <- "wget ftp://prism.nacse.org/monthly/tmean/"
# mid_ppt <- "/PRISM_ppt_stable_4kmM3_"
# mid_tmp <- "/PRISM_tmean_stable_4kmM2_"
# end <- "_all_bil.zip"
# 
# cat(c("#!/usr/bin/env bash", paste0(base_ppt, year, mid_ppt, year, end)),
#     file = "bash/fetch_prism_ppt.sh", sep = "\n")
# cat(c("#!/usr/bin/env bash", paste0(base_tmp, year, mid_tmp, year, end)),
#     file = "bash/fetch_prism_tmp.sh", sep = "\n")



# function to get clim data from raster file for given set of coordinates
fetch_prism <- function(file_tmp, file_ppt, spp) {
  prism_tmp <- raster(x = file_tmp)
  prism_ppt <- raster(x = file_ppt)
  if (!is.data.frame(spp)) spp <- spp[[1]]
  coordinates(spp) <- c('Lon', 'Lat')
  pbase <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  proj4string(spp) <- pbase
  spp <- spTransform(spp, crs(prism_tmp))
  spp$ppt <- vapply(raster::extract(prism_ppt, spp, buffer = 4000), mean, na.rm = TRUE, numeric(1))
  spp$tmp <- vapply(raster::extract(prism_tmp, spp, buffer = 4000), mean, na.rm = TRUE, numeric(1))
  # spp$tmp <- prism_tmp[cellFromXY(prism_tmp, spp)]
  # spp$ppt <- prism_ppt[cellFromXY(prism_ppt, spp)]
  # spp$ppt[spp$ppt > 65000] <- NA_real_
  return(as_tibble(spp))
}


# prism raster files
prism_files <- paste0("prism/", list.files("prism"))
month_files <- prism_files[grepl("[[:digit:]]{6}", prism_files)]
bil_files <- month_files[grepl(".bil$", month_files)]
files_ppt <- bil_files[grepl("ppt", bil_files)]
files_tmp <- bil_files[grepl("tmean", bil_files)]


# coordinates for each species/pop of interest
# spp_df <- read_csv("data/clim/species_coords.csv")

spp_df <- tibble(
  SpeciesAuthor = "Astragalus_scaphoides_2",
  MatrixPopulation = c("Haynes Creek", "McDevitt Creek", "Sheep Corral Gulch"),
  Lat = c(45.005, 44.923, 45.105),
  Lon = c(-113.702, -113.729 , -113.045)
)

spp_df <- tibble(SpeciesAuthor = "Cirsium_undulatum",
                 MatrixPopulation = "Hays",
                 Lat = 38.8, Lon = -99.3)


# get climate data from all raster files for all species of interest
df_clim <- tibble(file_tmp = files_tmp, file_ppt = files_ppt) %>% 
  mutate(Date = map_chr(files_tmp, ~ strsplit(.x, "_")[[1]][5])) %>% 
  mutate(Year = as.integer(substr(Date, 1, 4))) %>% 
  mutate(Month = as.integer(substr(Date, 5,6))) %>% 
  group_by(Year, Month) %>% 
  do(fetch_prism(.$file_tmp, .$file_ppt, spp_df)) %>% 
  ungroup() %>% 
  arrange(SpeciesAuthor, MatrixPopulation, Year, Month) %>% 
  filter(!(is.na(tmp) & is.na(ppt)))

write.csv(df_clim, "data/clim/species_clim_prism_hays_buffer.csv", row.names = FALSE)


