# S04: generate study-level sampling distributions for analysis 1.

# libraries ----
source("code/setup.R")
setup_packages(c("tidyverse", "Rcompadre", "popbio"))
source("code/functions.R")
set.seed(5654)


# load compadre data ----
compadre <- load_compadre(corrected = TRUE)
compadre <- compadre %>%
  mutate(
    MatrixStartYear = suppressWarnings(as.integer(MatrixStartYear)),
    MatrixEndYear = suppressWarnings(as.integer(MatrixEndYear)),
    YearPublication = suppressWarnings(as.integer(YearPublication))
  )

if (!("Observation" %in% names(compadre)) && ("Observations" %in% names(compadre))) {
  compadre <- compadre %>% mutate(Observation = as.character(Observations))
}
if (!("Observation" %in% names(compadre))) {
  compadre <- compadre %>% mutate(Observation = NA_character_)
}


# Load data from Ellis et al. (2012) ----
ellis_data <- read.table("data/raw/ellis_2012/Transition_Matrices.txt",
  sep = "\t",
  header = TRUE, stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(matA = map(Mx, string_to_mat)) %>%
  mutate(matU = map(Tmx, string_to_mat)) %>%
  mutate(matF = map2(matA, matU, ~ .x - .y)) %>%
  mutate(N = map(Nx, nx_to_vec))


# Draw from MPM sampling distributions by study ----

# Aschero ----
spp <- "Prosopis_ﬂexuosa"
aschero_n <- read_csv("data/derived/studies/aschero_n.csv")

aschero <- compadre %>%
  filter(SpeciesAuthor == spp, MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  mutate(N = list(aschero_n$N))

sd_aschero <- aschero %>%
  mutate(simU = map2(matU, N, ~ sim_U_wrapper(matU = .x, N = .y, nsim = 1000))) %>%
  mutate(simF = map2(matF, N, ~ sim_F_wrapper(matF = .x, N = .y, nsim = 1000))) %>%
  as_tibble() %>%
  select(MatrixPopulation, simU, simF)

aschero_out <- compadre %>%
  filter(SpeciesAuthor == spp, MatrixTreatment == "Unmanipulated") %>%
  left_join(sd_aschero)

save(aschero_out, file = "data/derived/analysis_cache/sd_aschero.RData")

dataf <- aschero %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = FALSE)


# Kiviniemi ----
spp <- "Agrimonia_eupatoria"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

kiviniemi_n <- read_csv("data/derived/studies/kiviniemi_n.csv") %>%
  group_by(SpeciesAccepted, MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

kiviniemi <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(kiviniemi_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

dataf <- kiviniemi %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Summarize boundary survival estimates and small stage sample sizes ----
sd_files <- Sys.glob("data/derived/analysis_cache/sd_*.RData")

all_sd <- map_dfr(sd_files, function(path) {
  out <- rdata_load2(path)
  as_tibble(out) %>%
    mutate(source_file = basename(path))
})

population_base <- all_sd %>%
  mutate(
    MatrixDimension = suppressWarnings(as.integer(MatrixDimension)),
    surv = map(mat, ~ colSums(Rcompadre::matU(.x))),
    stage_n = map(N, as.numeric)
  ) %>%
  transmute(
    source_file,
    Authors,
    YearPublication,
    SpeciesAccepted,
    MatrixPopulation,
    OrganismType,
    MatrixDimension,
    surv,
    stage_n
  )

stage_diag <- population_base %>%
  mutate(stage_tbl = map2(surv, stage_n, function(s, n) {
    if (length(n) == 0 || length(n) != length(s)) {
      n <- rep(NA_real_, length(s))
    }

    tibble(
      stage_index = seq_along(s),
      survival = as.numeric(s),
      n_stage = as.numeric(n)
    )
  })) %>%
  select(-surv, -stage_n) %>%
  unnest(stage_tbl) %>%
  mutate(
    boundary_zero = survival == 0,
    boundary_one = survival == 1,
    boundary_any = boundary_zero | boundary_one,
    small_n20 = n_stage < 20,
    small_n50 = n_stage < 50
  )

population_diag <- stage_diag %>%
  group_by(
    Authors,
    YearPublication,
    SpeciesAccepted,
    MatrixPopulation,
    OrganismType,
    MatrixDimension
  ) %>%
  summarize(
    n_stages = dplyr::n(),
    n_stages_with_N = sum(!is.na(n_stage)),
    any_boundary = any(boundary_any, na.rm = TRUE),
    any_boundary_one = any(boundary_one, na.rm = TRUE),
    prop_boundary = mean(boundary_any, na.rm = TRUE),
    prop_boundary_one = mean(boundary_one, na.rm = TRUE),
    any_small_n20 = any(small_n20, na.rm = TRUE),
    any_small_n50 = any(small_n50, na.rm = TRUE),
    prop_small_n20 = mean(small_n20, na.rm = TRUE),
    prop_small_n50 = mean(small_n50, na.rm = TRUE),
    min_stage_n = suppressWarnings(min(n_stage, na.rm = TRUE)),
    median_stage_n = median(n_stage, na.rm = TRUE),
    mean_stage_n = mean(n_stage, na.rm = TRUE),
    mean_survival = mean(survival, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(
      c(min_stage_n, median_stage_n, mean_stage_n),
      ~ ifelse(is.infinite(.x), NA_real_, .x)
    ),
    dim_group = if_else(MatrixDimension <= 5, "3-5", "6+")
  )

population_with_counts <- population_diag %>%
  filter(n_stages_with_N > 0)

by_life_form <- population_with_counts %>%
  group_by(OrganismType) %>%
  summarize(
    group_type = "Life form",
    n_populations = dplyr::n(),
    pct_any_boundary = mean(any_boundary) * 100,
    pct_any_boundary_one = mean(any_boundary_one) * 100,
    pct_any_small_n20 = mean(any_small_n20) * 100,
    pct_any_small_n50 = mean(any_small_n50) * 100,
    mean_prop_boundary = mean(prop_boundary),
    mean_prop_small_n20 = mean(prop_small_n20),
    mean_min_stage_n = mean(min_stage_n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(group_value = OrganismType)

by_dimension_group <- population_with_counts %>%
  group_by(dim_group) %>%
  summarize(
    group_type = "Matrix dimension",
    n_populations = dplyr::n(),
    pct_any_boundary = mean(any_boundary) * 100,
    pct_any_boundary_one = mean(any_boundary_one) * 100,
    pct_any_small_n20 = mean(any_small_n20) * 100,
    pct_any_small_n50 = mean(any_small_n50) * 100,
    mean_prop_boundary = mean(prop_boundary),
    mean_prop_small_n20 = mean(prop_small_n20),
    mean_min_stage_n = mean(min_stage_n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(group_value = dim_group)

by_dimension <- population_with_counts %>%
  group_by(MatrixDimension) %>%
  summarize(
    n_populations = dplyr::n(),
    pct_any_boundary = mean(any_boundary) * 100,
    pct_any_boundary_one = mean(any_boundary_one) * 100,
    pct_any_small_n20 = mean(any_small_n20) * 100,
    pct_any_small_n50 = mean(any_small_n50) * 100,
    mean_prop_boundary = mean(prop_boundary),
    mean_prop_small_n20 = mean(prop_small_n20),
    mean_min_stage_n = mean(min_stage_n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(MatrixDimension)

overall <- tibble(
  n_populations_all = nrow(population_diag),
  n_populations_with_counts = nrow(population_with_counts),
  n_stage_estimates = nrow(stage_diag),
  pct_pop_with_any_boundary = mean(population_with_counts$any_boundary) * 100,
  pct_pop_with_any_boundary_one = mean(population_with_counts$any_boundary_one) * 100,
  pct_pop_with_any_small_n20 = mean(population_with_counts$any_small_n20) * 100,
  pct_pop_with_any_small_n50 = mean(population_with_counts$any_small_n50) * 100,
  corr_prop_boundary_prop_small_n20 = cor(
    population_with_counts$prop_boundary,
    population_with_counts$prop_small_n20
  ),
  corr_matrix_dimension_prop_boundary = cor(
    population_with_counts$MatrixDimension,
    population_with_counts$prop_boundary
  ),
  corr_matrix_dimension_prop_small_n20 = cor(
    population_with_counts$MatrixDimension,
    population_with_counts$prop_small_n20
  )
)

group_summary <- bind_rows(
  by_life_form,
  by_dimension_group
)

write_csv(
  stage_diag,
  "data/derived/analysis_cache/boundary_smallN_stage_summary.csv"
)
write_csv(
  population_diag,
  "data/derived/analysis_cache/boundary_smallN_population_summary.csv"
)
write_csv(
  group_summary,
  "data/derived/analysis_cache/boundary_smallN_group_summary.csv"
)
write_csv(
  by_dimension,
  "data/derived/analysis_cache/boundary_smallN_by_dimension.csv"
)
write_csv(
  overall,
  "data/derived/analysis_cache/boundary_smallN_overall.csv"
)


npool <- kiviniemi %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

kiviniemi_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  filter(MatrixPopulation %in% c("A", "B")) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(kiviniemi_out, file = "data/derived/analysis_cache/sd_kiviniemi.RData")


# Satterthwaite ----
spp <- "Eriogonum_longifolium_var._gnaphalifolium_2"
pop <- "Unburned"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

satterthwaite_n <- read_csv("data/derived/studies/satterthwaite_n.csv") %>%
  mutate(Nf = N) %>%
  mutate(Nu = ifelse(Pool, 0, N)) %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(Nu = list(Nu), Nf = list(Nf)) %>%
  ungroup()

satterthwaite <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(satterthwaite_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

has_counts <- function(x) {
  !is.null(x) && length(x) > 0 && !all(is.na(x))
}

satterthwaite_tbl <- satterthwaite %>%
  as_tibble() %>%
  mutate(
    has_Nu = map_lgl(Nu, has_counts),
    has_Nf = map_lgl(Nf, has_counts)
  )

satterthwaite_excluded <- satterthwaite_tbl %>%
  filter(!has_Nu | !has_Nf) %>%
  distinct(MatrixPopulation, MatrixStartYear, MatrixEndYear, Authors, DOI_ISBN)

satterthwaite_latest_summary <- satterthwaite_tbl %>%
  group_by(MatrixPopulation) %>%
  summarize(
    latest_n_mats = n(),
    latest_year_min = suppressWarnings(min(as.integer(MatrixStartYear), na.rm = TRUE)),
    latest_year_max = suppressWarnings(max(as.integer(MatrixStartYear), na.rm = TRUE)),
    missing_count_rows = sum(!has_Nu | !has_Nf),
    .groups = "drop"
  )

old_compadre_candidates <- list.files(
  "data/raw/compadre",
  pattern = "^COMPADRE_v\\.X\\.X\\.X_pre_case1_check_.*\\.RData$",
  full.names = TRUE
)
old_compadre_path <- if (length(old_compadre_candidates) > 0) {
  old_compadre_candidates[which.max(file.info(old_compadre_candidates)$mtime)]
} else {
  NA_character_
}

if (!is.na(old_compadre_path) && file.exists(old_compadre_path)) {
  old_compadre <- cdb_fetch(old_compadre_path)
  old_satterthwaite_summary <- old_compadre %>%
    filter(SpeciesAuthor == spp) %>%
    filter(MatrixComposite == "Individual") %>%
    filter(MatrixTreatment == "Unmanipulated") %>%
    cdb_unnest() %>%
    as_tibble() %>%
    group_by(MatrixPopulation) %>%
    summarize(
      old_n_mats = n(),
      old_year_min = suppressWarnings(min(as.integer(MatrixStartYear), na.rm = TRUE)),
      old_year_max = suppressWarnings(max(as.integer(MatrixStartYear), na.rm = TRUE)),
      .groups = "drop"
    )
} else {
  old_satterthwaite_summary <- tibble(
    MatrixPopulation = character(),
    old_n_mats = integer(),
    old_year_min = integer(),
    old_year_max = integer()
  )
}

satterthwaite_accounting <- full_join(
  old_satterthwaite_summary,
  satterthwaite_latest_summary,
  by = "MatrixPopulation"
) %>%
  mutate(
    excluded_in_latest = MatrixPopulation %in% satterthwaite_excluded$MatrixPopulation,
    exclusion_reason = ifelse(
      excluded_in_latest,
      "No matching stage-count rows in data/derived/studies/satterthwaite_n.csv",
      NA_character_
    )
  ) %>%
  arrange(MatrixPopulation)

write_csv(
  satterthwaite_accounting,
  "data/derived/studies/satterthwaite_population_accounting.csv"
)

write_csv(
  satterthwaite_excluded,
  "data/derived/studies/satterthwaite_excluded_latest.csv"
)

npool <- satterthwaite_tbl %>%
  filter(has_Nu) %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(Nu)), .groups = "drop")

satterthwaite_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  filter(MatrixPopulation %in% npool$MatrixPopulation) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(satterthwaite_out, file = "data/derived/analysis_cache/sd_satterthwaite.RData")

dataf <- satterthwaite %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)

# Andrello ----
spp <- "Eryngium_alpinum"
pop <- "PRD" # DES, BER, BOU, PRA, PRB, PRC, PRD


andrello_n <- read_csv("data/derived/studies/andrello_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

andrello <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  cdb_flag("check_zero_U") %>%
  filter(check_zero_U == FALSE) %>%
  left_join(andrello_n, by = c("MatrixPopulation", "MatrixStartYear")) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- andrello %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

andrello_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  population_matrices_from_available(preferred = c("Mean", "Individual")) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(andrello_out, file = "data/derived/analysis_cache/sd_andrello.RData")

dataf <- andrello %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Liatris_scariosa ----
spp <- "Liatris_scariosa"
# Ellis: LISC_0, LISC_1, LISC_2
# Comp: "Lisc 0", "Lisc 1", "Lisc 2"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse()

lisc_n <- ellis_data %>%
  filter(SPP == "LISC") %>%
  mutate(MatrixPopulation = case_when(
    POP == "LISC_0" ~ "Lisc 0",
    POP == "LISC_1" ~ "Lisc 1",
    POP == "LISC_2" ~ "Lisc 2"
  )) %>%
  select(MatrixPopulation, MatrixStartYear = YR, N)

lisc <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(lisc_n, by = c("MatrixPopulation", "MatrixStartYear")) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- lisc %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

lisc_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(!grepl(";", MatrixPopulation)) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(lisc_out, file = "data/derived/analysis_cache/sd_lisc.RData")

dataf <- lisc %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Cirsium_pitcheri_4 ----
# Compadre has CiPi 1, CiPi 2, CiPi 3; Ellis has CIPI_1, CIPI_2, CIPI_3, CIPI_4
# I think Cirsium_pitcheri_6 from Bell et al 2013, corresponds to CIPI 4 from
#  Ellis et al 2012 (Cirsium_pitcheri_4), but they use diff stage classes
spp <- "Cirsium_pitcheri_4"
pop <- "CIPI_3" # Ellis: "CIPI_1", "CIPI_2", "CIPI_3"
pop_comp <- "CiPi 3" # Comp: "CiPi 1", "CiPi 2", "CiPi 3"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse()

cipi_n <- ellis_data %>%
  filter(SPP == "CIPI") %>%
  mutate(MatrixPopulation = case_when(
    POP == "CIPI_1" ~ "CiPi 1",
    POP == "CIPI_2" ~ "CiPi 2",
    POP == "CIPI_3" ~ "CiPi 3"
  )) %>%
  mutate(PU = map(N, function(x) ifelse(x == 0, TRUE, FALSE))) %>%
  mutate(PF = map(N, function(x) ifelse(x == 0, TRUE, FALSE))) %>%
  select(MatrixPopulation, MatrixStartYear = YR, N, PU, PF)

cipi <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(cipi_n, by = c("MatrixPopulation", "MatrixStartYear")) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- cipi %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

cipi_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(!grepl(";", MatrixPopulation)) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(cipi_out, file = "data/derived/analysis_cache/sd_cipi.RData")

dataf <- cipi %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Scanga ----
spp <- "Trollius_laxus_2"
pop <- c("CfCh", "Cb", "EEFF", "H66cont", "MM", "T")

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixPopulation %in% pop) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

scanga_n <- read_csv("data/derived/studies/scanga_n.csv") %>%
  rename(MatrixPopulation = Group) %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(N)) %>%
  ungroup()

scanga <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixPopulation %in% pop) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(scanga_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- scanga %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

scanga_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixPopulation %in% pop) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(scanga_out, file = "data/derived/analysis_cache/sd_scanga.RData")

dataf <- scanga %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Lazaro ----
spp <- "Dioon_merolae"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

lazaro_n <- read_csv("data/derived/studies/lazaro_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

lazaro <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(lazaro_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- lazaro %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

lazaro_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  population_matrices_from_available(preferred = c("Mean", "Individual")) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(lazaro_out, file = "data/derived/analysis_cache/sd_lazaro.RData")

dataf <- lazaro %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Arroyo ----
spp <- "Neobuxbaumia_polylopha"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

arroyo_n <- read_csv("data/derived/studies/arroyo_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

arroyo <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(arroyo_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- arroyo %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

arroyo_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(arroyo_out, file = "data/derived/analysis_cache/sd_arroyo.RData")

dataf <- arroyo %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Plank ----
spp <- "Trillium_persistens"
# "Battle Creek", "Moccasin Creek", "Moody Creek", "Panther Creek"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

plank_n <- read_csv("data/derived/studies/plank_n.csv") %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(N)) %>%
  ungroup()

plank <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  filter(!grepl("fecundity", MatrixPopulation, ignore.case = TRUE)) %>%
  cdb_unnest() %>%
  left_join(plank_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- plank %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)), .groups = "drop")

# moody and panther had 0 seedlings... use pooled value of 5 instead
npool <- npool %>%
  mutate(
    N = case_when(
      MatrixPopulation %in% c("Moody Creek", "Panther Creek") ~ map(N, ~ {
        .x[1] <- 5
        .x
      }),
      TRUE ~ N
    )
  )

plank_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  filter(!grepl("fecundity", MatrixPopulation, ignore.case = TRUE)) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(plank_out, file = "data/derived/analysis_cache/sd_plank.RData")

dataf <- plank %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Jolls ----
spp <- "Cirsium_pitcheri_8"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

jolls_n <- read_csv("data/derived/studies/jolls_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

jolls <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(jolls_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup() %>%
  mutate(MatrixPopulation = ifelse(MatrixStartYear <= 2000, "1995", "2005"))

npool <- jolls %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

jolls_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  mutate(MatrixPopulation = ifelse(MatrixStartYear <= 2000, "1995", "2005")) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(jolls_out, file = "data/derived/analysis_cache/sd_jolls.RData")

dataf <- jolls %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Torres ----
spp <- "Agave_potatorum"
# "Xochiltepec", "Machiche"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

torres_n <- read_csv("data/derived/studies/torres_n.csv") %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(N)) %>%
  ungroup()

torres <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(torres_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- torres %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

torres_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(torres_out, file = "data/derived/analysis_cache/sd_torres.RData")

dataf <- torres %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Andrieu ----
spp <- "Paeonia_officinalis"
pops <- c("Open habitat", "Woodland")
# "Open habitat", "Woodland", (Managed habitat doesn't have pooled)

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixPopulation %in% pops) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

andrieu_n <- read_csv("data/derived/studies/andrieu_n.csv") %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(N)) %>%
  ungroup()

andrieu <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixPopulation %in% pops) %>%
  filter(MatrixComposite == "Pooled") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(andrieu_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- andrieu %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

andrieu_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixPopulation %in% pops) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  population_matrices_from_available(preferred = c("Pooled", "Mean", "Individual")) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(andrieu_out, file = "data/derived/analysis_cache/sd_andrieu.RData")

dataf <- andrieu %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Eriksson ----
spp <- "Plantago_media"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

eriksson_n <- read_csv("data/derived/studies/eriksson_n.csv") %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(N)) %>%
  ungroup()

eriksson <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(eriksson_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- eriksson %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

eriksson_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  slice(-grep(";", MatrixPopulation)) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(eriksson_out, file = "data/derived/analysis_cache/sd_eriksson.RData")

dataf <- eriksson %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Astragalus_scaphoides_2, Haynes Creek, Sheep Corral Gulch, McDevitt Creek ----
# sometimes 0 fecund
# negative relationship between fecundity and sample size
spp <- "Astragalus_scaphoides_2"
# ASSC_haynes, ASSC_sheep, ASSC_mcdevi
# Haynes Creek, Sheep Corral Gulch, McDevitt Creek

compadre %>%
  filter(SpeciesAuthor == "Astragalus_scaphoides_2") %>%
  cdb_glimpse("MatrixComposite")

assc_n <- ellis_data %>%
  filter(SPP == "ASSC") %>%
  mutate(MatrixPopulation = case_when(
    POP == "ASSC_haynes" ~ "Haynes Creek",
    POP == "ASSC_sheep" ~ "Sheep Corral Gulch",
    POP == "ASSC_mcdevi" ~ "McDevitt Creek"
  )) %>%
  mutate(PU = map(N, function(x) ifelse(x == 0, TRUE, FALSE))) %>%
  mutate(PF = map(N, function(x) ifelse(x == 0, TRUE, FALSE))) %>%
  select(MatrixPopulation, MatrixStartYear = YR, N, PU, PF)

assc <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(assc_n, by = c("MatrixPopulation", "MatrixStartYear")) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- assc %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

assc_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  population_matrices_from_available(preferred = c("Mean", "Individual")) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(assc_out, file = "data/derived/analysis_cache/sd_assc.RData")

dataf <- assc %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Lemke ----
spp <- "Trollius_europaeus"
# "HAS; JAG", "RDGm; GTH; SPW; NEV", "RDGab; JAGab"
# *NOTE* "HAS; JAG" has N = 0 repro, so only use for surv analyses ----

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  filter(MatrixComposite == "Pooled") %>%
  filter(Observation == "Pooled by habitat and year") %>%
  cdb_glimpse(c("Observation", "MatrixComposite"))

lemke_n <- read_csv("data/derived/studies/lemke_n.csv") %>%
  group_by(MatrixPopulation, Observation, MatrixStartYear) %>%
  summarize(N = list(N))

lemke <- compadre %>%
  filter(
    SpeciesAuthor == spp,
    MatrixComposite == "Pooled",
    Observation == "Pooled by habitat and year"
  ) %>%
  mutate(Observation = as.character(Observation)) %>%
  cdb_unnest() %>%
  left_join(lemke_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- lemke %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

lemke_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  filter(MatrixComposite == "Pooled") %>%
  filter(Observation == "Pooled by habitat and year") %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(lemke_out, file = "data/derived/analysis_cache/sd_lemke.RData")

dataf <- lemke %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Toledo ----
# matrix values based on bootstrapping, so won't necessarily match
spp <- "Tillandsia_butzii"
# "San Antonio, Veracruz"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse()

toledo_n <- read_csv("data/derived/studies/toledo_n.csv") %>%
  group_by(MatrixStartYear) %>%
  summarize(N = list(N))

toledo_raw <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated")

if (nrow(as_tibble(toledo_raw)) == 0) {
  old_compadre_candidates <- list.files(
    "data/raw/compadre",
    pattern = "^COMPADRE_v\\.X\\.X\\.X_pre_case1_check_.*\\.RData$",
    full.names = TRUE
  )
  old_compadre_path <- if (length(old_compadre_candidates) > 0) {
    old_compadre_candidates[which.max(file.info(old_compadre_candidates)$mtime)]
  } else {
    NA_character_
  }

  old_n_individual <- NA_integer_
  if (!is.na(old_compadre_path) && file.exists(old_compadre_path)) {
    old_compadre <- cdb_fetch(old_compadre_path)
    old_n_individual <- old_compadre %>%
      filter(SpeciesAuthor == spp) %>%
      filter(MatrixComposite == "Individual") %>%
      filter(MatrixTreatment == "Unmanipulated") %>%
      as_tibble() %>%
      nrow()
  }

  toledo_exclusion <- tibble(
    study = "Toledo",
    SpeciesAuthor = spp,
    reason = "No matching individual unmanipulated matrices in active COMPADRE version",
    old_n_individual = old_n_individual,
    new_n_individual = 0L
  )
  write_csv(
    toledo_exclusion,
    "data/derived/studies/toledo_excluded_latest.csv"
  )

  warning(
    "Skipping Toledo block: no individual unmanipulated Tillandsia_butzii matrices in active COMPADRE."
  )
} else {

toledo <- toledo_raw %>%
  cdb_unnest() %>%
  left_join(toledo_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- toledo %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

toledo_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(toledo_out, file = "data/derived/analysis_cache/sd_toledo.RData")

dataf <- toledo %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)
}


# Crone ----
spp <- "Balsamorhiza_sagittata"
# "Mount Jumbo"
# fecundity based on number of flowers, not plants

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse()

crone_n <- read_csv("data/derived/studies/crone_n.csv") %>%
  group_by(MatrixStartYear) %>%
  summarize(N = list(N))

crone <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(crone_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- crone %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

crone_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(crone_out, file = "data/derived/analysis_cache/sd_crone.RData")

dataf <- crone %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Dostalek ----
spp <- "Dracocephalum_austriacum_2"
# Cisarska rokle (C1), Haknovec (C2), Kodska stena (C3)
# Zadielsky kamen (S1), Domicke skrapy (S2), Zelezne vrata (S3)
# seed survival constant across sites (N = 97)

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse()

dostalek_n <- read_csv("data/derived/studies/dostalek_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

dostalek <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(dostalek_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- dostalek %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

dostalek_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  population_matrices_from_available(preferred = c("Mean", "Individual")) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(dostalek_out, file = "data/derived/analysis_cache/sd_dostalek.RData")

dataf <- dostalek %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Evju ----
spp <- "Viola_biflora"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse()

evju_n <- read_csv("data/derived/studies/evju_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

evju <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(evju_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- evju %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

evju_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(evju_out, file = "data/derived/analysis_cache/sd_evju.RData")

dataf <- evju %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Flores ----
spp <- "Mammillaria_huitzilopochtli"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse()

flores_n <- read_csv("data/derived/studies/flores_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

flores <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(flores_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- flores %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

flores_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(flores_out, file = "data/derived/analysis_cache/sd_flores.RData")

dataf <- flores %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Shryock ----
spp <- "Pediocactus_bradyi"

shryock_n <- read_csv("data/derived/studies/shryock_n.csv") %>%
  mutate(N = pmap(list(S1, S2, S3), ~ c(..1, ..2, ..3))) %>%
  mutate(PU = map(N, function(x) ifelse(x == 0, TRUE, FALSE))) %>%
  select(MatrixPopulation, MatrixStartYear, N, PU)

shryock <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  cdb_unnest() %>%
  left_join(shryock_n, by = c("MatrixPopulation", "MatrixStartYear")) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  mutate(posF = list(mat_mean(matF) > 0)) %>%
  ungroup()

npool <- shryock %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

shryock_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  population_matrices_from_available(preferred = c("Mean", "Individual")) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(shryock_out, file = "data/derived/analysis_cache/sd_shryock.RData")

dataf <- shryock %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Csergo ----
spp <- "Saponaria_bellidifolia"

csergo_n <- read_csv("data/derived/studies/csergo_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup() %>%
  select(MatrixPopulation, MatrixStartYear, N)

csergo <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  cdb_unnest() %>%
  left_join(csergo_n, by = c("MatrixPopulation", "MatrixStartYear")) %>%
  group_by(MatrixPopulation) %>%
  mutate(
    posU = list(mat_mean(matU) > 0),
    posF = list(mat_mean(matF) > 0)
  ) %>%
  ungroup()

npool <- csergo %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

csergo_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(!grepl(";", MatrixPopulation)) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(csergo_out, file = "data/derived/analysis_cache/sd_csergo.RData")

dataf <- csergo %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Raghu ----
spp <- "Lantana_camara_2"

raghu_n <- read_csv("data/derived/studies/raghu_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup() %>%
  select(MatrixPopulation, MatrixStartYear, N)

raghu <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(
    MatrixComposite == "Individual",
    MatrixTreatment == "Unmanipulated"
  ) %>%
  cdb_unnest() %>%
  left_join(raghu_n, by = c("MatrixPopulation", "MatrixStartYear")) %>%
  group_by(MatrixPopulation) %>%
  mutate(
    posU = list(mat_mean(matU) > 0),
    posF = list(mat_mean(matF) > 0)
  ) %>%
  ungroup()

npool <- raghu %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

raghu_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean", MatrixTreatment == "Unmanipulated") %>%
  filter(!grepl(";", MatrixPopulation)) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(raghu_out, file = "data/derived/analysis_cache/sd_raghu.RData")

dataf <- raghu %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Martin ----
spp <- "Astragalus_peckii"

martin_n <- read_csv("data/derived/studies/martin_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup() %>%
  select(MatrixPopulation, MatrixStartYear, N)

martin <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(
    MatrixComposite == "Individual",
    MatrixTreatment == "Unmanipulated"
  ) %>%
  cdb_unnest() %>%
  left_join(martin_n, by = c("MatrixPopulation", "MatrixStartYear")) %>%
  group_by(MatrixPopulation) %>%
  mutate(
    posU = list(mat_mean(matU) > 0),
    posF = list(mat_mean(matF) > 0)
  ) %>%
  ungroup()

npool <- martin %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

martin_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  filter(MatrixComposite == "Mean") %>%
  filter(!grepl(";", MatrixPopulation)) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(martin_out, file = "data/derived/analysis_cache/sd_martin.RData")

dataf <- martin %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Law ----
# single pooled value of fecundity, from 38 individs
spp <- "Saussurea_medusa"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse()

law_n <- read_csv("data/derived/studies/law_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

law <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(law_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  ungroup()

npool <- law %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

npool$N[[1]][5] <- 38

law_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(law_out, file = "data/derived/analysis_cache/sd_law.RData")

dataf <- law %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)

# Jacquemyns ----
spp <- "Orchis_purpurea"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse()

jacq_n <- read_csv("data/derived/studies/jacquemyns_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

jacq <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(jacq_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(
    posU = list(mat_mean(matU) > 0),
    posF = list(mat_mean(matF) > 0)
  ) %>%
  ungroup()

npool <- jacq %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

jacq_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(!grepl(";", MatrixPopulation)) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(jacq_out, file = "data/derived/analysis_cache/sd_jacq.RData")

dataf <- jacq %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Portela ----
spp <- "Astrocaryum_aculeatissimum"

portela_n <- read_csv("data/derived/studies/portela_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup() %>%
  select(MatrixPopulation, MatrixStartYear, N)

portela <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(
    MatrixComposite == "Individual",
    MatrixTreatment == "Unmanipulated"
  ) %>%
  cdb_unnest() %>%
  left_join(portela_n, by = c("MatrixPopulation", "MatrixStartYear")) %>%
  group_by(MatrixPopulation) %>%
  mutate(
    posU = list(mat_mean(matU) > 0),
    posF = list(mat_mean(matF) > 0)
  ) %>%
  ungroup()

npool <- portela %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

portela_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean", MatrixTreatment == "Unmanipulated") %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(portela_out, file = "data/derived/analysis_cache/sd_portela.RData")

dataf <- portela %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Lopez-mata ----
spp <- "Pinus_maximartinezii"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse()

lopez_n <- read_csv("data/derived/studies/lopez_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

lopez <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(lopez_n) %>%
  mutate(posU = map(matU, ~ .x > 0))

npool <- lopez %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

lopez_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(lopez_out, file = "data/derived/analysis_cache/sd_lopez.RData")

dataf <- lopez %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Auestad ----
spp <- "Pimpinella_saxifraga"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse()

auestad_n <- read_csv("data/derived/studies/auestad_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

auestad <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(auestad_n) %>%
  mutate(posU = map(matU, ~ .x > 0)) %>%
  mutate(posF = map(matF, ~ .x > 0))

npool <- auestad %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

auestad_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(auestad_out, file = "data/derived/analysis_cache/sd_auestad.RData")

dataf <- auestad %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)


# Dias Segura ----
spp <- "Lophophora_diffusa"

compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse()

dias_n <- read_csv("data/derived/studies/dias_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>%
  summarize(N = list(N)) %>%
  ungroup()

dias <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Individual") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_unnest() %>%
  left_join(dias_n) %>%
  group_by(MatrixPopulation) %>%
  mutate(posU = list(mat_mean(matU) > 0)) %>%
  ungroup()

npool <- dias %>%
  as_tibble() %>%
  group_by(MatrixPopulation) %>%
  summarize(N = list(pool_counts(N)))

dias_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(!grepl(";", MatrixPopulation)) %>%
  left_join(npool) %>%
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>%
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) 

save(dias_out, file = "data/derived/analysis_cache/sd_dias.RData")

dataf <- dias %>%
  cdb_metadata() %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)
dataf <- unique(dataf)
mdata <- paste(dataf$Authors, dataf$YearPublication, dataf$Journal, dataf$DOI_ISBN, dataf$SpeciesAccepted, sep = ", ")

write(mdata, file = "data/derived/studies/_data_sources.csv", append = TRUE)
