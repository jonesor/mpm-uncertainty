

### libraries
library(tidyverse)
library(Rcompadre)
library(Rage)
library(popbio)
library(gridExtra)
source("R/functions.R")


### load compadre data
compadre <- cdb_fetch("data/COMPADRE_v.X.X.X_Corrected.RData")


### Load data from Ellis et al. (2012)
ellis_data <- read.table('data/ellis/Transition_Matrices.txt', sep = '\t',
                         header = TRUE, stringsAsFactors = FALSE) %>%
  as_tibble() %>% 
  mutate(matA = lapply(Mx, stringToMat)) %>% 
  mutate(matU = lapply(Tmx, stringToMat)) %>% 
  mutate(matF = mapply(function(a, b) a - b, matA, matU, SIMPLIFY = F)) %>% 
  mutate(N = lapply(Nx, NxToVec))





### Draw from MPM sampling distributions by study ##############################

### Aschero
spp <- "Prosopis_ﬂexuosa"

aschero_n <- read_csv("data/studies/aschero_n.csv")

aschero <- compadre %>% 
  filter(SpeciesAuthor == spp, MatrixTreatment == "Unmanipulated") %>% 
  cdb_unnest() %>% 
  mutate(N = list(aschero_n$N))

sd_aschero <- aschero %>% 
  mutate(simU = map2(matU, N, ~ SimMatUWrapper(matU = .x, N = .y, nsim = 1000))) %>% 
  mutate(simF = map2(matF, N, ~ SimMatFWrapper(matF = .x, N = .y, nsim = 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, simU, simF)

aschero_out <- compadre %>% 
  filter(SpeciesAuthor == spp, MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_aschero)

# save(aschero_out, file = "analysis/sd_aschero.RData")



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

# sampling distribution
sd_kiviniemi <- kiviniemi %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatFWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

kiviniemi_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  filter(MatrixPopulation %in% c("A", "B")) %>% 
  left_join(sd_kiviniemi)

# save(kiviniemi_out, file = "analysis/sd_kiviniemi.RData")



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

# sampling distribution
sd_satterthwaite <- satterthwaite %>% 
  mutate(simU = pmap(list(matU, posU, Nu), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, Nf), ~ SimMatFWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(popbio::mean.list(simU, na.rm = TRUE)),
            simF = list(popbio::mean.list(simF, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

satterthwaite_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_satterthwaite)

# save(satterthwaite_out, file = "analysis/sd_satterthwaite.RData")



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

# sampling distribution
sd_andrello <- andrello %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatFWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

andrello_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  slice(-grepl(";", MatrixPopulation)) %>% 
  left_join(sd_andrello)

# save(andrello_out, file = "analysis/sd_andrello.RData")



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

# sampling distribution
sd_lisc <- lisc %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatFWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(popbio::mean.list(simU, na.rm = TRUE)),
            simF = list(popbio::mean.list(simF, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

lisc_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(!grepl(";", MatrixPopulation)) %>% 
  left_join(sd_lisc)

# save(lisc_out, file = "analysis/sd_lisc.RData")



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

# sampling distribution
sd_cipi <- cipi %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatFWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(popbio::mean.list(simU, na.rm = TRUE)),
            simF = list(popbio::mean.list(simF, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

cipi_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(!grepl(";", MatrixPopulation)) %>% 
  left_join(sd_cipi)

# save(cipi_out, file = "analysis/sd_cipi.RData")



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

# sampling distribution
sd_scanga <- scanga %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

scanga_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixPopulation %in% pop) %>% 
  filter(MatrixComposite == "Individual") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_scanga)

# save(scanga_out, file = "analysis/sd_scanga.RData")



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

# sampling distribution
sd_lazaro <- lazaro %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

lazaro_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  slice(-grep(";", MatrixPopulation)) %>% 
  left_join(sd_lazaro)

# save(lazaro_out, file = "analysis/sd_lazaro.RData")



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

# sampling distribution
sd_arroyo <- arroyo %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatFWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

arroyo_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_arroyo)

# save(arroyo_out, file = "analysis/sd_arroyo.RData")



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

# moody and panther had 0 seedlings... use pooled value of 5 instead
plank$N[[3]][1] <- 5
plank$N[[4]][1] <- 5

# sampling distribution
sd_plank <- plank %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean2(simU)),
            simF = list(mat_mean2(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

plank_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Individual") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_plank)

# save(plank_out, file = "analysis/sd_plank.RData")



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
  mutate(MatrixPopulation = ifelse(MatrixStartYear <= 2000, 1995, 2005))

# sampling distribution
sd_jolls <- jolls %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatFWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

jolls_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_jolls, by = c("MatrixStartYear" = "MatrixPopulation"))

# save(jolls_out, file = "analysis/sd_jolls.RData")



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

# sampling distribution
sd_torres <- torres %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatFWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

torres_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Individual") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_torres)

# save(torres_out, file = "analysis/sd_torres.RData")



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

andrieu_n <- read_csv("data/studies/andrieu_n_pooled.csv") %>% 
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

# sampling distribution
sd_andrieu <- andrieu %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

andrieu_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixPopulation %in% pops) %>%
  filter(MatrixComposite == "Pooled") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_andrieu)

# save(andrieu_out, file = "analysis/sd_andrieu.RData")



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

# sampling distribution
sd_eriksson <- eriksson %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

eriksson_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  slice(-grep(";", MatrixPopulation)) %>% 
  left_join(sd_eriksson)

# save(eriksson_out, file = "analysis/sd_eriksson.RData")



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
  left_join(assc_n, by = c("MatrixPopulation", "MatrixStartYear"))

mat_dim <- ncol(assc$matA[[1]])

# fit bayesian models, and save posterior samples
trans_sim <- assc %>% 
  as_tibble() %>% 
  group_by(MatrixPopulation) %>% 
  do(bind_rows(lapply(seq_len(mat_dim),
                      tr_bayes_wrap,
                      matU = .$matU, matF = .$matF, N = .$N,
                      poolU = .$PU, poolF = .$PF))) %>% 
  ungroup()

# df of all transitions, possible and not
cast_dim <- function(dim) {
  dim <- unique(dim[[1]])
  expand.grid(col = seq_len(dim),
              row = seq_len(dim),
              rep = seq_len(1000),
              stringsAsFactors = FALSE)
}

tr_df <- assc %>%
  as_tibble() %>%
  select(MatrixPopulation, year = MatrixStartYear) %>%
  group_by(MatrixPopulation) %>%
  mutate(year = as.integer(as.factor(year))) %>%
  ungroup() %>%
  group_by(MatrixPopulation, year) %>%
  do(cast_dim(mat_dim)) %>%
  ungroup()

trans_sim_full <- trans_sim %>%
  select(MatrixPopulation, year) %>%
  unique() %>%
  mutate(SpeciesAuthor = spp) %>%
  full_join(tr_df) %>%
  select(-SpeciesAuthor) %>%
  left_join(trans_sim) %>%
  arrange(MatrixPopulation, year, rep, col, row)

sd_assc <- trans_sim_full %>%
  mutate(tr_U = ifelse(is.na(tr_U), 0, tr_U),
         tr_F = ifelse(is.na(tr_F), 0, tr_F)) %>%
  group_by(MatrixPopulation, rep, col, row) %>%
  summarize(tr_U = mean(tr_U),
            tr_F = mean(tr_F)) %>%
  ungroup() %>%
  group_by(MatrixPopulation, rep) %>%
  summarize(simU = list(MakeMat(tr_U, mat_dim)),
            simF = list(MakeMat(tr_F, mat_dim))) %>%
  ungroup() %>%
  group_by(MatrixPopulation) %>%
  summarize(simU = list(simU),
            simF = list(simF)) %>%
  mutate(SpeciesAuthor = spp)

assc_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>%
  filter(MatrixComposite == "Mean") %>%
  filter(MatrixTreatment == "Unmanipulated") %>%
  slice(-grep(";", MatrixPopulation)) %>% 
  left_join(sd_assc)

# save(assc_out, file = "analysis/sd_assc.RData")



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

# sampling distribution
sd_lemke <- lemke %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean2(simU)),
            simF = list(mat_mean2(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

lemke_out <- compadre %>%
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  filter(MatrixComposite == "Pooled") %>% 
  filter(Observation == "Pooled by habitat and year") %>%
  left_join(sd_lemke)

# save(lemke_out, file = "analysis/sd_lemke.RData")



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
  ungroup()

# sampling distribution
sd_toledo <- toledo %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

toledo_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_toledo)

# save(toledo_out, file = "analysis/sd_toledo.RData")



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
  ungroup()

# sampling distribution
sd_crone <- crone %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

crone_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_crone)

# save(crone_out, file = "analysis/sd_crone.RData")



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
  ungroup()

# sampling distribution
sd_dostalek <- dostalek %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

dostalek_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_dostalek)

# save(dostalek_out, file = "analysis/sd_dostalek.RData")



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
  ungroup()

# sampling distribution
sd_evju <- evju %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

evju_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_evju)

# save(evju_out, file = "analysis/sd_evju.RData")



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
  ungroup()

# sampling distribution
sd_flores <- flores %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU)),
            simF = list(mat_mean(simF))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

flores_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_flores)

# save(flores_out, file = "analysis/sd_flores.RData")



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
  ungroup()

# sampling distribution
sd_shryock <- shryock %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(popbio::mean.list(simU, na.rm = TRUE)),
            simF = list(popbio::mean.list(simF, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

shryock_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(!grepl(";", MatrixPopulation)) %>% 
  left_join(sd_shryock)

# save(shryock_out, file = "analysis/sd_shryock.RData")



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

# sampling distribution
sd_csergo <- csergo %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(popbio::mean.list(simU, na.rm = TRUE)),
            simF = list(popbio::mean.list(simF, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

csergo_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(!grepl(";", MatrixPopulation)) %>% 
  left_join(sd_csergo)

# save(csergo_out, file = "analysis/sd_csergo.RData")



### Keller
spp <- "Leontopodium_alpinum"

keller_n <- read_csv("data/studies/keller_n.csv") %>%
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  summarize(N = list(N)) %>% 
  ungroup() %>% 
  select(MatrixPopulation, MatrixStartYear, N)

keller <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Individual") %>% 
  cdb_unnest() %>% 
  mutate(matF = map2(matF, matC, ~ .x + .y)) %>% 
  mutate(matC = map(matC, ~ matrix(0, nrow(.x), ncol(.x)))) %>% 
  left_join(keller_n, by = c("MatrixPopulation", "MatrixStartYear")) %>% 
  group_by(MatrixPopulation) %>% 
  mutate(posU = list(mat_mean(matU) > 0),
         posF = list(mat_mean(matF) > 0)) %>% 
  ungroup()

# sampling distribution
sd_keller <- keller %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatFWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(popbio::mean.list(simU, na.rm = TRUE)),
            simF = list(popbio::mean.list(simF, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

keller_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(!grepl(";", MatrixPopulation)) %>% 
  left_join(sd_keller)

for (i in 1:nrow(keller_out)) {
  keller_out$mat[[i]]@matF <- keller_out$mat[[i]]@matF + keller_out$mat[[i]]@matC
  keller_out$mat[[i]]@matC[keller_out$mat[[i]]@matC > 0] <- 0
}

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

# sampling distribution
sd_raghu <- raghu %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatFWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(popbio::mean.list(simU)),
            simF = list(popbio::mean.list(simF, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

raghu_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean", MatrixTreatment == "Unmanipulated") %>% 
  filter(!grepl(";", MatrixPopulation)) %>% 
  left_join(sd_raghu)

# save(raghu_out, file = "analysis/sd_raghu.RData")



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
  mutate(posU = list(mat_mean(matU) > 0)) %>% 
  ungroup()

# sampling distribution
sd_martin <- martin %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(popbio::mean.list(simU, na.rm = TRUE)),
            simF = list(popbio::mean.list(simF, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

martin_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(!grepl(";", MatrixPopulation)) %>% 
  left_join(sd_martin)

# save(martin_out, file = "analysis/sd_martin.RData")



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

# sampling distribution
sd_law <- law %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(mat_mean(simU))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU))

law_fec <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  cdb_unnest() %>% 
  mutate(N = list(c(0, 0, 0, 0, 74))) %>% 
  mutate(posF = list(mat_mean(matF) > 0)) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatFWrapper(..1, ..2, ..3, 1000)))

sd_law$simF <- law_fec$simF

law_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_law)

# save(law_out, file = "analysis/sd_law.RData")



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
  mutate(posU = list(mat_mean(matU) > 0)) %>% 
  ungroup()

# sampling distribution
sd_jacq <- jacq %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(popbio::mean.list(simU, na.rm = TRUE)),
            simF = list(popbio::mean.list(simF, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

jacq_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(!grepl(";", MatrixPopulation)) %>% 
  left_join(sd_jacq)

# save(jacq_out, file = "analysis/sd_jacq.RData")



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

# sampling distribution
sd_portela <- portela %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = pmap(list(matF, posF, N), ~ SimMatFWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(popbio::mean.list(simU)),
            simF = list(popbio::mean.list(simF, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

portela_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean", MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_portela)

# save(portela_out, file = "analysis/sd_portela.RData")



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

# sampling distribution
sd_lopez <- lopez %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

lopez_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  left_join(sd_lopez)

# save(lopez_out, file = "analysis/sd_lopez.RData")



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
  mutate(posU = map(matU, ~ .x > 0))

# sampling distribution
sd_auestad <- auestad %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU) %>% 
  unnest() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU))

auestad_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(sd_auestad) %>% 
  mutate(simF = map(matF(mat), ~ replicate(1000, .x, simplify = FALSE)))

# save(auestad_out, file = "analysis/sd_auestad.RData")



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

# sampling distribution
sd_dias <- dias %>% 
  mutate(simU = pmap(list(matU, posU, N), ~ SimMatUWrapper(..1, ..2, ..3, 1000))) %>% 
  mutate(simF = map(matF, ~ replicate(1000, .x, simplify = FALSE))) %>% 
  as_tibble() %>% 
  select(MatrixPopulation, MatrixStartYear, simU, simF) %>% 
  unnest() %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation, rep) %>% 
  summarize(simU = list(popbio::mean.list(simU, na.rm = TRUE)),
            simF = list(popbio::mean.list(simF, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(MatrixPopulation) %>% 
  summarize(simU = list(simU),
            simF = list(simF))

dias_out <- compadre %>% 
  filter(SpeciesAuthor == spp) %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(!grepl(";", MatrixPopulation)) %>% 
  left_join(sd_dias)

# save(dias_out, file = "analysis/sd_dias.RData")

