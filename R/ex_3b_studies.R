

### libraries
library(tidyverse)
library(Rcompadre)
library(popbio)
source("R/functions.R")


### load compadre data
compadre <- cdb_fetch("data/COMPADRE_v.X.X.X_Corrected.RData")


### Load data from Ellis et al. (2012)
ellis_data <- read.table("data/ellis/Transition_Matrices.txt", sep = "\t",
                         header = TRUE, stringsAsFactors = FALSE) %>%
  as_tibble() %>% 
  mutate(matA = lapply(Mx, string_to_mat)) %>% 
  mutate(matU = lapply(Tmx, string_to_mat)) %>% 
  mutate(matF = mapply(function(a, b) a - b, matA, matU, SIMPLIFY = FALSE)) %>% 
  mutate(N = lapply(Nx, nx_to_vec))




### Draw from MPM sampling distributions by study ##############################

### Aschero
spp <- "Prosopis_ﬂexuosa"
aschero_n <- read_csv("data/studies/aschero_n.csv")

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

save(aschero_out, file = "analysis/sd_aschero.RData")



### Kiviniemi
spp <- "Agrimonia_eupatoria"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == 'Unmanipulated') %>%
  cdb_glimpse("MatrixComposite")

kiviniemi_n <- read_csv("data/studies/kiviniemi_n.csv") %>% 
  group_by(SpeciesAccepted, MatrixPopulation, MatrixStartYear) %>% 
  summarize(N = list(N)) %>% 
  ungroup()

kiviniemi <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == 'Individual') %>% 
  filter(MatrixTreatment == 'Unmanipulated') %>% 
  cdb_unnest() %>% 
  left_join(kiviniemi_n) %>% 
  group_by(MatrixPopulation) %>% 
  mutate(posU = list(mat_mean(matU) > 0)) %>% 
  mutate(posF = list(mat_mean(matF) > 0)) %>% 
  ungroup()

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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(kiviniemi_out, file = "analysis/sd_kiviniemi.RData")



### Satterthwaite
spp <- "Eriogonum_longifolium_var._gnaphalifolium_2"
pop <- "Unburned"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == 'Unmanipulated') %>%
  cdb_glimpse("MatrixComposite")

satterthwaite_n <- read_csv("data/studies/satterthwaite_n.csv") %>% 
  mutate(Nf = N) %>% 
  mutate(Nu = ifelse(Pool, 0, N)) %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  summarize(Nu = list(Nu), Nf = list(Nf)) %>% 
  ungroup()

satterthwaite <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == 'Individual') %>% 
  filter(MatrixTreatment == 'Unmanipulated') %>% 
  cdb_unnest() %>% 
  left_join(satterthwaite_n) %>% 
  group_by(MatrixPopulation) %>% 
  mutate(posU = list(mat_mean(matU) > 0)) %>% 
  mutate(posF = list(mat_mean(matF) > 0)) %>% 
  ungroup()

npool <- satterthwaite %>% 
  as_tibble() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(N = list(pool_counts(Nu)))

satterthwaite_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(npool) %>% 
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(satterthwaite_out, file = "analysis/sd_satterthwaite.RData")



## Andrello
spp <- "Eryngium_alpinum"
pop <- "PRD" # DES, BER, BOU, PRA, PRB, PRC, PRD

# load("data/studies/raw/andrello_matrices.RData")
# 
# andrello_n <- apply(M, c(2, 3, 4), sum) %>%
#   as.data.frame() %>% 
#   as_tibble() %>% 
#   rownames_to_column("Stage") %>% 
#   gather(group, N, -Stage) %>% 
#   mutate(MatrixStartYear = as.integer(substr(group, 5, 9))) %>% 
#   mutate(MatrixPopulation = substr(group, 1, 3)) %>% 
#   mutate(SpeciesAuthor = spp) %>% 
#   select(SpeciesAuthor, Stage, MatrixPopulation, MatrixStartYear, N)
# 
# write.csv(andrello_n, "data/studies/andrello_n.csv", row.names = FALSE)

andrello_n <- read_csv("data/studies/andrello_n.csv") %>% 
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
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  slice(-grepl(";", MatrixPopulation)) %>% 
  left_join(npool) %>% 
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(andrello_out, file = "analysis/sd_andrello.RData")



### Liatris_scariosa
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
    POP == "LISC_2" ~ "Lisc 2")) %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(lisc_out, file = "analysis/sd_lisc.RData")



### Cirsium_pitcheri_4
# Compadre has CiPi 1, CiPi 2, CiPi 3; Ellis has CIPI_1, CIPI_2, CIPI_3, CIPI_4
# I think Cirsium_pitcheri_6 from Bell et al 2013, corresponds to CIPI 4 from
#  Ellis et al 2012 (Cirsium_pitcheri_4), but they use diff stage classes
spp <- "Cirsium_pitcheri_4"
pop <- "CIPI_3" # Ellis: "CIPI_1", "CIPI_2", "CIPI_3"
pop_comp <- "CiPi 3" # Comp: "CiPi 1", "CiPi 2", "CiPi 3"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == 'Mean') %>% 
  filter(MatrixTreatment == 'Unmanipulated') %>%
  cdb_glimpse()

cipi_n <- ellis_data %>%
  filter(SPP == 'CIPI') %>% 
  mutate(MatrixPopulation = case_when(
    POP == "CIPI_1" ~ "CiPi 1",
    POP == "CIPI_2" ~ "CiPi 2",
    POP == "CIPI_3" ~ "CiPi 3")) %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(cipi_out, file = "analysis/sd_cipi.RData")



### Scanga
spp <- "Trollius_laxus_2"
pop <- c("CfCh", "Cb", "EEFF", "H66cont", "MM", "T")

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixPopulation %in% pop) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

scanga_n <- read_csv("data/studies/scanga_n.csv") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(scanga_out, file = "analysis/sd_scanga.RData")



### Lazaro
spp <- "Dioon_merolae"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse("MatrixComposite")

lazaro_n <- read_csv("data/studies/lazaro_n.csv") %>% 
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
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  slice(-grep(";", MatrixPopulation)) %>% 
  left_join(npool) %>% 
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(lazaro_out, file = "analysis/sd_lazaro.RData")



### Arroyo
spp <- "Neobuxbaumia_polylopha"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

arroyo_n <- read_csv("data/studies/arroyo_n.csv") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(arroyo_out, file = "analysis/sd_arroyo.RData")



### Plank
spp <- "Trillium_persistens"
# "Battle Creek", "Moccasin Creek", "Moody Creek", "Panther Creek"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse("MatrixComposite")

plank_n <- read_csv("data/studies/plank_n.csv") %>% 
  group_by(MatrixPopulation) %>% 
  summarize(N = list(N)) %>% 
  ungroup()

plank <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Individual") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_unnest() %>% 
  left_join(plank_n) %>% 
  group_by(MatrixPopulation) %>% 
  mutate(posU = list(mat_mean(matU) > 0)) %>% 
  mutate(posF = list(mat_mean(matF) > 0)) %>% 
  ungroup()

npool <- plank %>% 
  as_tibble() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(N = list(pool_counts(N)))

# moody and panther had 0 seedlings... use pooled value of 5 instead
npool$N[[3]][1] <- 5
npool$N[[4]][1] <- 5

plank_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Individual") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(npool) %>% 
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(plank_out, file = "analysis/sd_plank.RData")



### Jolls
spp <- "Cirsium_pitcheri_8"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

jolls_n <- read_csv("data/studies/jolls_n.csv") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(jolls_out, file = "analysis/sd_jolls.RData")



### Torres
spp <- "Agave_potatorum"
# "Xochiltepec", "Machiche"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse("MatrixComposite")

torres_n <- read_csv("data/studies/torres_n.csv") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(torres_out, file = "analysis/sd_torres.RData")



##### Andrieu
# some stages use pooled for every single year... figure out how to assess
spp <- "Paeonia_officinalis"
pops <- c("Open habitat", "Woodland")
# "Open habitat", "Woodland", (Managed habitat doesn't have pooled)

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixPopulation %in% pops) %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

andrieu_n <- read_csv("data/studies/andrieu_n.csv") %>% 
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
  filter(MatrixComposite == "Pooled") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(npool) %>% 
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(andrieu_out, file = "analysis/sd_andrieu.RData")



##### Eriksson
# some stages use pooled for every single year... figure out how to assess
spp <- "Plantago_media"
# pop <- "Site B" # Site A, Site B

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>%
  cdb_glimpse("MatrixComposite")

eriksson_n <- read_csv("data/studies/eriksson_n.csv") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(eriksson_out, file = "analysis/sd_eriksson.RData")



### Astragalus_scaphoides_2, Haynes Creek, Sheep Corral Gulch, McDevitt Creek
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
    POP == "ASSC_mcdevi" ~ "McDevitt Creek")) %>% 
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
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  filter(!grepl(";", MatrixPopulation)) %>% 
  left_join(npool) %>% 
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(assc_out, file = "analysis/sd_assc.RData")



#### Lemke
spp <- "Trollius_europaeus"
# "HAS; JAG", "RDGm; GTH; SPW; NEV", "RDGab; JAGab"
### *NOTE* "HAS; JAG" has N = 0 repro, so only use for surv analyses

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  filter(MatrixComposite == "Pooled") %>% 
  filter(Observation == "Pooled by habitat and year") %>%
  cdb_glimpse(c("Observation", "MatrixComposite"))

lemke_n <- read_csv("data/studies/lemke_n.csv") %>% 
  group_by(MatrixPopulation, Observation, MatrixStartYear) %>% 
  summarize(N = list(N))

lemke <- compadre %>% 
  filter(SpeciesAuthor == spp,
         MatrixComposite == "Pooled",
         Observation == "Pooled by habitat and year") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(lemke_out, file = "analysis/sd_lemke.RData")



### Toledo
# matrix values based on bootstrapping, so won't necessarily match
spp <- "Tillandsia_butzii"
# "San Antonio, Veracruz"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse()

toledo_n <- read_csv("data/studies/toledo_n.csv") %>% 
  group_by(MatrixStartYear) %>% 
  summarize(N = list(N))

toledo <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Individual") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(toledo_out, file = "analysis/sd_toledo.RData")



### Crone
spp <- "Balsamorhiza_sagittata"
# "Mount Jumbo"
# fecundity based on number of flowers, not plants

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse()

crone_n <- read_csv("data/studies/crone_n.csv") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(crone_out, file = "analysis/sd_crone.RData")



### Dostalek
spp <- "Dracocephalum_austriacum_2"
# Cisarska rokle (C1), Haknovec (C2), Kodska stena (C3)
# Zadielsky kamen (S1), Domicke skrapy (S2), Zelezne vrata (S3)
# seed survival constant across sites (N = 97)

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse()

dostalek_n <- read_csv("data/studies/dostalek_n.csv") %>% 
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
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(npool) %>% 
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(dostalek_out, file = "analysis/sd_dostalek.RData")



### Evju
spp <- "Viola_biflora"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse()

evju_n <- read_csv("data/studies/evju_n.csv") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(evju_out, file = "analysis/sd_evju.RData")



### Flores
spp <- "Mammillaria_huitzilopochtli"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse()

flores_n <- read_csv("data/studies/flores_n.csv") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(flores_out, file = "analysis/sd_flores.RData")



### Shryock
spp <- "Pediocactus_bradyi"

shryock_n <- read_csv("data/studies/shryock_n.csv") %>%
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
  filter(MatrixComposite == "Mean") %>% 
  filter(!grepl(";", MatrixPopulation)) %>% 
  left_join(npool) %>% 
  mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(shryock_out, file = "analysis/sd_shryock.RData")



### Csergo
spp <- "Saponaria_bellidifolia"

csergo_n <- read_csv("data/studies/csergo_n.csv") %>%
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
  mutate(posU = list(mat_mean(matU) > 0),
         posF = list(mat_mean(matF) > 0)) %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(csergo_out, file = "analysis/sd_csergo.RData")



# ### Keller
# spp <- "Leontopodium_alpinum"
# 
# keller_n <- read_csv("data/studies/keller_n.csv") %>%
#   group_by(MatrixPopulation, MatrixStartYear) %>% 
#   summarize(N = list(N)) %>% 
#   ungroup() %>% 
#   select(MatrixPopulation, MatrixStartYear, N)
# 
# keller <- compadre %>% 
#   filter(SpeciesAuthor == spp) %>% 
#   filter(MatrixComposite == "Individual") %>% 
#   cdb_unnest() %>% 
#   mutate(matF = map2(matF, matC, ~ .x + .y)) %>% 
#   mutate(matC = map(matC, ~ matrix(0, nrow(.x), ncol(.x)))) %>% 
#   left_join(keller_n, by = c("MatrixPopulation", "MatrixStartYear")) %>% 
#   group_by(MatrixPopulation) %>% 
#   mutate(posU = list(mat_mean(matU) > 0),
#          posF = list(mat_mean(matF) > 0)) %>% 
#   ungroup()
# 
# npool <- keller %>% 
#   as_tibble() %>% 
#   group_by(MatrixPopulation) %>% 
#   summarize(N = list(pool_counts(N)))
# 
# keller_out <- compadre %>% 
#   filter(SpeciesAuthor == spp) %>% 
#   filter(MatrixComposite == "Mean") %>% 
#   filter(!grepl(";", MatrixPopulation)) %>% 
#   left_join(npool) %>% 
#   mutate(simU = pmap(list(matU(mat), N), ~ sim_U_wrapper(..1, N = ..2, nsim = 1000))) %>% 
#   mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
#   select(-N)
# 
# for (i in 1:nrow(keller_out)) {
#   keller_out$mat[[i]]@matF <- keller_out$mat[[i]]@matF + keller_out$mat[[i]]@matC
#   keller_out$mat[[i]]@matC[keller_out$mat[[i]]@matC > 0] <- 0
# }
# 
# save(keller_out, file = "analysis/sd_keller.RData")



### Raghu
spp <- "Lantana_camara_2"

raghu_n <- read_csv("data/studies/raghu_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  summarize(N = list(N)) %>% 
  ungroup() %>% 
  select(MatrixPopulation, MatrixStartYear, N)

raghu <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Individual",
         MatrixTreatment == "Unmanipulated") %>% 
  cdb_unnest() %>% 
  left_join(raghu_n, by = c("MatrixPopulation", "MatrixStartYear")) %>% 
  group_by(MatrixPopulation) %>% 
  mutate(posU = list(mat_mean(matU) > 0),
         posF = list(mat_mean(matF) > 0)) %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(raghu_out, file = "analysis/sd_raghu.RData")



### Martin
spp <- "Astragalus_peckii"

martin_n <- read_csv("data/studies/martin_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  summarize(N = list(N)) %>% 
  ungroup() %>% 
  select(MatrixPopulation, MatrixStartYear, N)

martin <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Individual",
         MatrixTreatment == "Unmanipulated") %>% 
  cdb_unnest() %>% 
  left_join(martin_n, by = c("MatrixPopulation", "MatrixStartYear")) %>% 
  group_by(MatrixPopulation) %>% 
  mutate(posU = list(mat_mean(matU) > 0),
         posF = list(mat_mean(matF) > 0)) %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(martin_out, file = "analysis/sd_martin.RData")



### Law
# single pooled value of fecundity, from 38 individs
spp <- "Saussurea_medusa"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse()

law_n <- read_csv("data/studies/law_n.csv") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(law_out, file = "analysis/sd_law.RData")



### Jacquemyns
spp <- "Orchis_purpurea"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse()

jacq_n <- read_csv("data/studies/jacquemyns_n.csv") %>% 
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
  mutate(posU = list(mat_mean(matU) > 0),
         posF = list(mat_mean(matF) > 0)) %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(jacq_out, file = "analysis/sd_jacq.RData")



### Portela
spp <- "Astrocaryum_aculeatissimum"

portela_n <- read_csv("data/studies/portela_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  summarize(N = list(N)) %>% 
  ungroup() %>% 
  select(MatrixPopulation, MatrixStartYear, N)

portela <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Individual",
         MatrixTreatment == "Unmanipulated") %>% 
  cdb_unnest() %>% 
  left_join(portela_n, by = c("MatrixPopulation", "MatrixStartYear")) %>% 
  group_by(MatrixPopulation) %>% 
  mutate(posU = list(mat_mean(matU) > 0),
         posF = list(mat_mean(matF) > 0)) %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(portela_out, file = "analysis/sd_portela.RData")



### Lopez-mata
spp <- "Pinus_maximartinezii"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse()

lopez_n <- read_csv("data/studies/lopez_n.csv") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(lopez_out, file = "analysis/sd_lopez.RData")



### Auestad
spp <- "Pimpinella_saxifraga"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse()

auestad_n <- read_csv("data/studies/auestad_n.csv") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(auestad_out, file = "analysis/sd_auestad.RData")



### Dias Segura
spp <- "Lophophora_diffusa"

compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  cdb_glimpse()

dias_n <- read_csv("data/studies/dias_n.csv") %>% 
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
  mutate(simF = pmap(list(matF(mat), N), ~ sim_F_wrapper(..1, N = ..2, nsim = 1000))) %>% 
  select(-N)

save(dias_out, file = "analysis/sd_dias.RData")

