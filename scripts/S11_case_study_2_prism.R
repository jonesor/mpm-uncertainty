# S11: download/extract PRISM climate rasters and build climate inputs for analysis 2.

# libraries ----
source("code/setup.R")
setup_packages(c("tidyverse", "sf", "terra", "Rcompadre", "prism"))
source("code/functions.R")


# prism raster files ----
# download location is configured by setup_prism() in code/setup.R
prism_dir <- setup_prism("data/raw/prism")

maybe_download_prism <- function(years) {
  if (length(list.files(prism_dir, pattern = "\\.bil$", recursive = TRUE)) > 0) {
    return(invisible(NULL))
  }
  prism::get_prism_monthlys(type = "ppt", years = years, mon = 1:12, keepZip = FALSE)
  prism::get_prism_monthlys(type = "tmean", years = years, mon = 1:12, keepZip = FALSE)
}

bil_files <- list.files(prism_dir, pattern = "\\.bil$", recursive = TRUE, full.names = TRUE)
month_files <- bil_files[grepl("[[:digit:]]{6}", bil_files)]
files_ppt <- month_files[grepl("ppt", month_files)]
files_tmp <- month_files[grepl("tmean", month_files)]


# coordinates for each species/pop of interest ----
coords_path <- "data/derived/climate/species_coords.csv"
if (!file.exists(coords_path)) {
  compadre <- load_compadre(corrected = TRUE)
  comp_sub <- compadre %>%
    filter(
      MatrixComposite == "Individual",
      MatrixTreatment == "Unmanipulated",
      ProjectionInterval == "1",
      MatrixCaptivity == "W"
    )
  comp_time_series <- comp_sub %>%
    as_tibble() %>%
    filter(!is.na(Lon) & !is.na(Lat)) %>%
    group_by(SpeciesAuthor) %>%
    mutate(n_year = length(unique(MatrixStartYear))) %>%
    ungroup() %>%
    filter(n_year >= 5) %>%
    group_by(SpeciesAuthor, MatrixPopulation) %>%
    summarize(
      Lon = unique(Lon)[1],
      Lat = unique(Lat)[1],
      n_year = unique(n_year),
      .groups = "drop"
    )
  write_csv(comp_time_series, coords_path)
}

# target sites for climate extraction ----
# Keep this list explicit so analysis sites are reproducible and easy to extend.
target_sites <- tibble(
  SpeciesAuthor = c(
    "Silene_spaldingii",
    "Astragalus_scaphoides_2",
    "Astragalus_scaphoides_2",
    "Astragalus_scaphoides_2"
  ),
  MatrixPopulation = c(
    "Eureka",
    "Haynes Creek",
    "Sheep Corral Gulch",
    "McDevitt Creek"
  )
)

spp_df <- read_csv(coords_path) %>%
  semi_join(target_sites, by = c("SpeciesAuthor", "MatrixPopulation"))

if (length(bil_files) == 0) {
  compadre <- load_compadre(corrected = TRUE)
  years <- compadre %>%
    as_tibble() %>%
    semi_join(target_sites, by = c("SpeciesAuthor", "MatrixPopulation")) %>%
    filter(!is.na(MatrixStartYear)) %>%
    pull(MatrixStartYear) %>%
    unique() %>%
    sort()
  maybe_download_prism(years)
  bil_files <- list.files(prism_dir, pattern = "\\.bil$", recursive = TRUE, full.names = TRUE)
  month_files <- bil_files[grepl("[[:digit:]]{6}", bil_files)]
  files_ppt <- month_files[grepl("ppt", month_files)]
  files_tmp <- month_files[grepl("tmean", month_files)]
}


# get climate data from all raster files for all species of interest ----
df_ppt <- tibble(file_ppt = files_ppt) %>%
  mutate(Date = map_chr(file_ppt, ~ strsplit(basename(.x), "_")[[1]][5]))
df_tmp <- tibble(file_tmp = files_tmp) %>%
  mutate(Date = map_chr(file_tmp, ~ strsplit(basename(.x), "_")[[1]][5]))

df_clim <- inner_join(df_tmp, df_ppt, by = "Date") %>%
  mutate(Year = as.integer(substr(Date, 1, 4))) %>%
  mutate(Month = as.integer(substr(Date, 5, 6))) %>%
  group_by(Year, Month) %>%
  do(fetch_prism(.$file_tmp, .$file_ppt, spp_df)) %>%
  ungroup() %>%
  arrange(SpeciesAuthor, MatrixPopulation, Year, Month) %>%
  filter(!(is.na(tmp) & is.na(ppt)))

write.csv(df_clim, "data/derived/climate/species_clim_prism.csv", row.names = FALSE)
