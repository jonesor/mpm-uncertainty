
# libraries
library(rgdal)
library(raster)
library(readr)
library(tibble)
library(dplyr)
library(purrr)


## create shell scripts to download chelsa climate rasters
# year <- 1979:2013
# month <- formatC(1:12, width = 2, flag = "0")
# ym <- expand.grid(month = month, year = year)
# 
# tmp <- "wget https://www.wsl.ch/lud/chelsa/data/timeseries/tmean/CHELSA_tmean"
# ppt <- "wget https://www.wsl.ch/lud/chelsa/data/timeseries/prec/CHELSA_prec"
# end <- "V1.2.1.tif"
# 
# cat(c("#!/usr/bin/env bash", paste(tmp, ym$year, ym$month, end, sep = "_")),
#     file = "fetch_chelsa_tmp.sh", sep = "\n")
# cat(c("#!/usr/bin/env bash", paste(ppt, ym$year, ym$month, end, sep = "_")),
#     file = "fetch_chelsa_ppt.sh", sep = "\n")


# function to get clim data from raster file for given set of coordinates
fetch_chelsa <- function(file_tmp, file_ppt, spp) {
  chelsa_tmp <- raster(x = file_tmp)
  chelsa_ppt <- raster(x = file_ppt)
  if (!is.data.frame(spp)) spp <- spp[[1]]
  coordinates(spp) <- c('Lon', 'Lat')
  pbase <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  proj4string(spp) <- pbase
  spp <- spTransform(spp, crs(chelsa_tmp))
  spp$tmp <- chelsa_tmp[cellFromXY(chelsa_tmp, spp)] / 10 - 273.15
  spp$ppt <- chelsa_ppt[cellFromXY(chelsa_ppt, spp)]
  # spp$ppt[spp$ppt > 65000] <- NA_real_
  return(as_tibble(spp))
}


# chelsa raster files
chelsa_files <- paste0("chelsa/", list.files("chelsa"))
files_tmp <- chelsa_files[grepl("tmean", chelsa_files)]
files_ppt <- chelsa_files[grepl("prec", chelsa_files)]

# coordinates for each species/pop of interest
spp_df <- read_csv("data/clim/compadre_coords.csv")

# get climate data from all raster files for all species of interest
df_clim <- tibble(file_tmp = files_tmp, file_ppt = files_ppt,
                  spp = list(spp_df)) %>% 
  mutate(Year = map_int(files_tmp, ~ as.integer(strsplit(.x, "_")[[1]][3]))) %>% 
  mutate(Month = map_int(files_tmp, ~ as.integer(strsplit(.x, "_")[[1]][4]))) %>% 
  group_by(Year, Month) %>% 
  do(fetch_chelsa(.$file_tmp, .$file_ppt, .$spp)) %>% 
  ungroup() %>% 
  arrange(SpeciesAuthor, MatrixPopulation, Year, Month)

write.csv(df_clim, "data/clim/compadre_chelsa.csv", row.names = FALSE)


