
### libraries
library(tidyverse)
library(Rcompadre)
library(Rage)
library(popbio)
library(gridExtra)
source("R/functions.R")



### read chelsa climate data
spp_clim <- read_csv("data/clim/species_clim_prism_hays.csv")

astrag_wx <- read_csv("data/clim/weather_astragalus.csv") %>% 
  select(Year = years, Month = months, MatrixPopulation = Site, astr_ppt = ppt,
         astr_tmp = tmean) %>% 
  mutate(MatrixPopulation = case_when(
    MatrixPopulation == "Haynes" ~ "Haynes Creek",
    MatrixPopulation == "McDevitt" ~ "McDevitt Creek",
    MatrixPopulation == "Sheep" ~ "Sheep Corral Gulch"))

out <- spp_clim %>% 
  filter(SpeciesAuthor == "Astragalus_scaphoides_2") %>% 
  left_join(astrag_wx) %>% 
  filter(!is.na(astr_ppt))





adler_ppt <- read_csv("data/clim/adler_monthly_ppt.csv") %>% 
  filter(YEAR >= 1960) %>%
  setNames(c("Year", 1:12)) %>% 
  gather(Month, ad_ppt, -Year) %>% 
  mutate(Month = as.integer(Month)) %>% 
  arrange(Year, Month) %>% 
  mutate(n = 1:n())

adler_clim <- read_csv("data/clim/adler_monthly_temp.csv") %>% 
  filter(YEAR >= 1960) %>%
  setNames(c("Year", 1:12)) %>% 
  gather(Month, ad_tmp, -Year) %>% 
  mutate(Month = as.integer(Month)) %>% 
  arrange(Year, Month) %>% 
  left_join(adler_ppt)

out <- spp_clim %>% 
  filter(SpeciesAuthor == "Cirsium_undulatum") %>% 
  right_join(adler_clim) %>% 
  filter(!is.na(tmp))

df_lim <- out %>% 
  group_by(Month) %>% 
  summarize(min_tmp = min(c(tmp, ad_tmp), na.rm = TRUE),
            min_ppt = min(c(ppt, ad_ppt), na.rm = TRUE),
            max_tmp = max(c(tmp, ad_tmp), na.rm = TRUE),
            max_ppt = max(c(ppt, ad_ppt), na.rm = TRUE))



ggplot(out, aes(tmp, ad_tmp)) +
  geom_blank(data = df_lim, aes(x = min_tmp, y = min_tmp)) +
  geom_blank(data = df_lim, aes(x = max_tmp, y = max_tmp)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.6) +
  geom_point(aes(color = Year)) +
  facet_wrap(~ Month, scales = "free") +
  theme(panel.grid = element_blank())

ggplot(out, aes(ppt, ad_ppt)) +
  geom_blank(data = df_lim, aes(x = min_ppt, y = min_ppt)) +
  geom_blank(data = df_lim, aes(x = max_ppt, y = max_ppt)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.6) +
  geom_point(aes(color = Year)) +
  facet_wrap(~ Month, scales = "free") +
  theme(panel.grid = element_blank())




# McDevitt Creek, Haynes Creek, Sheep Corral Gulch
pop <- "Haynes Creek"

df_lim <- out %>% 
  filter(MatrixPopulation == pop) %>% 
  group_by(Month) %>% 
  summarize(min_tmp = min(c(tmp, astr_tmp), na.rm = TRUE),
            min_ppt = min(c(ppt, astr_ppt), na.rm = TRUE),
            max_tmp = max(c(tmp, astr_tmp), na.rm = TRUE),
            max_ppt = max(c(ppt, astr_ppt), na.rm = TRUE))

ggplot(filter(out, MatrixPopulation == pop), aes(tmp, astr_tmp)) +
  geom_blank(data = df_lim, aes(x = min_tmp, y = min_tmp)) +
  geom_blank(data = df_lim, aes(x = max_tmp, y = max_tmp)) +
  geom_point(aes(color = Year)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.6) +
  facet_wrap(~ Month, scales = "free") +
  theme(panel.grid = element_blank())

ggplot(filter(out, MatrixPopulation == pop), aes(ppt, astr_ppt)) +
  geom_blank(data = df_lim, aes(x = min_ppt, y = min_ppt)) +
  geom_blank(data = df_lim, aes(x = max_ppt, y = max_ppt)) +
  geom_point(aes(color = Year)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.6) +
  facet_wrap(~ Month, scales = "free") +
  theme(panel.grid = element_blank())


out %>% 
  group_by(MatrixPopulation, Month) %>% 
  summarize(r2 = cor(ppt, ad_ppt, use = "complete.obs")^2) %>% 
  ungroup() %>% 
  as.data.frame()







