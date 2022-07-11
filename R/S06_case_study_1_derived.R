
### libraries
library(tidyverse)
library(popbio)
library(popdemo)
library(Rcompadre)
library(Rage)
source("R/functions.R")


### load compadre data
compadre <- cdb_fetch("data/COMPADRE_v.X.X.X_Corrected.RData")


### load study-specific sampling distribution files
sd_files <- paste0("analysis/", list.files("analysis"))
sd_files <- sd_files[grep("/sd_", sd_files)]


### bind sampling distributions into single tibble
mpm_draws <- cdb_bind_rows(lapply(sd_files, rdata_load)) %>% 
  mutate(id = as.factor(1:n())) %>% 
  cdb_unnest() %>% 
  mutate(any_repro = map_lgl(matF, ~ any(.x > 0))) %>% 
  filter(any_repro == TRUE) %>% # make sure some repro
  mutate(matU = map(matU, scale_U)) %>% 
  mutate(matA = pmap(list(matU, matF, matC), ~ ..1 + ..2 + ..3)) %>% 
  mutate(start = map_int(mat, Rcompadre::mpm_first_active)) %>% 
  mutate(rep_stages = map(matF, ~ colSums(.x) > 0)) %>% 
  mutate(exclude_stages = map(MatrixClassOrganized,
                              ~ ifelse(.x == "active", FALSE, TRUE)))

# summary counts
length(unique(mpm_draws$Authors))
length(unique(mpm_draws$MatrixPopulation))
length(unique(mpm_draws$SpeciesAccepted))

# n_pops by species
mpm_draws %>% 
  as_tibble() %>% 
  count(SpeciesAuthor) %>% 
  arrange(desc(n))


### point estimates for parameters of interest
pt_shape <- mpm_draws %>% 
  mutate(rep_prop1 = pmap(list(matU, start, rep_stages), Rage::mature_distrib)) %>% 
  mutate(lx4 = map2_dbl(matU, rep_prop1,
                        ~ Rage::mpm_to_lx(.x, .y, lx_crit = -1, xmax = 3)[4])) %>% 
  filter(lx4 > 0) %>% # check perennial
  mutate(q = map2_int(matU, rep_prop1, ~ Rage::qsd_converge(.x, .y, conv = 0.01, N = 1e5))) %>% 
  filter(q >= 3) %>%   # make sure at least 3 time steps
  mutate(lx = pmap(list(matU, rep_prop1, q),
                   ~ Rage::mpm_to_lx(..1, ..2, xmax = ..3), lx_crit = -1)) %>% 
  mutate(L_pt = map2_dbl(matU, rep_prop1, Rage::life_expect)) %>% 
  mutate(lx_min = map_dbl(lx, min)) %>% 
  mutate(S_pt = map_dbl(lx, Rage::shape_surv)) %>% 
  as_tibble() %>% 
  mutate(id_L = fct_reorder(fct_drop(id), L_pt)) %>%
  mutate(id_S = fct_reorder(fct_drop(id), S_pt))

pt_other <- mpm_draws %>% 
  mutate(loglam_pt = map_dbl(matA, ~ log(popbio::lambda(.x)))) %>% 
  mutate(damp_pt = map_dbl(matA, popbio::damping.ratio)) %>% 
  mutate(gen_pt = map2_dbl(matU, matF, Rage::gen_time)) %>%
  mutate(pmature_pt = pmap_dbl(list(matU, matF, start), Rage::mature_prob)) %>% 
  mutate(growth_pt = pmap_dbl(list(matU, exclude_stages),
                              ~ Rage::vr_growth(..1, exclude = ..2))) %>% 
  mutate(elast_pt = pmap_dbl(list(matU, matF, exclude_stages),
                             ~ perturb_cust(..1, ..2, exclude = ..3,
                                            type = "elasticity")$progr)) %>% 
  as_tibble() %>% 
  mutate(id_loglam = fct_reorder(fct_drop(id), loglam_pt)) %>%
  mutate(id_damp = fct_reorder(fct_drop(id), damp_pt)) %>%
  mutate(id_gen = fct_reorder(fct_drop(id), gen_pt)) %>%
  mutate(id_pmature = fct_reorder(fct_drop(id), pmature_pt)) %>%
  mutate(id_growth = fct_reorder(fct_drop(id), growth_pt)) %>%
  mutate(id_elast = fct_reorder(fct_drop(id), elast_pt))



### sampling distributions for derived parameters
sd_shape <- pt_shape %>%
  select(id, SpeciesAuthor, MatrixPopulation, simU, simF, q) %>%
  unnest(cols = c("simU", "simF")) %>%
  left_join(select(pt_shape, id, start, rep_stages), by = "id") %>%
  mutate(rep_prop1 = pmap(list(simU, start, rep_stages), Rage::mature_distrib)) %>%
  mutate(lx = pmap(list(simU, rep_prop1, q),
                   ~ Rage::mpm_to_lx(..1, ..2, xmax = ..3), lx_crit = -1)) %>%
  mutate(L = map2_dbl(simU, rep_prop1, Rage::life_expect)) %>%
  mutate(S = map_dbl(lx, Rage::shape_surv))

sd_other <- pt_other %>%
  select(id, SpeciesAuthor, MatrixPopulation, simU, simF) %>%
  unnest(cols = c("simU", "simF")) %>%
  left_join(select(pt_other, id, exclude_stages), by = "id") %>%
  mutate(simA = pmap(list(simU, simF), ~ ..1 + ..2)) %>%
  left_join(select(as_tibble(pt_other), id, start, rep_stages), by = "id") %>%
  mutate(loglam = map_dbl(simA, ~ log(popbio::lambda(.x)))) %>%
  mutate(damp = map_dbl(simA, popbio::damping.ratio)) %>%
  mutate(gen = map2_dbl(simU, simF, Rage::gen_time)) %>%
  mutate(pmature = pmap_dbl(list(simU, simF, start), Rage::mature_prob)) %>%
  mutate(growth = pmap_dbl(list(simU, exclude_stages),
                           ~ Rage::vr_growth(..1, exclude = ..2))) %>%
  mutate(elast = pmap_dbl(list(simU, simF, exclude_stages),
                          ~ perturb_cust(..1, ..2, exclude = ..3,
                                         type = "elasticity")$progr))


### write to file
pt_shape_out <- pt_shape %>%
  dplyr::select(id, SpeciesAuthor, MatrixPopulation, q,
                ends_with("_pt"), starts_with("id_"))

pt_other_out <- pt_other %>%
  dplyr::select(id, SpeciesAuthor, MatrixPopulation,
                ends_with("_pt"), starts_with("id_"))

sd_shape_out <- sd_shape %>%
  left_join(dplyr::select(pt_shape_out, id, starts_with("id_"))) %>%
  dplyr::select(id, SpeciesAuthor, MatrixPopulation, start, L, S, starts_with("id_"))

sd_other_out <- sd_other %>%
  left_join(dplyr::select(pt_other_out, id, starts_with("id_"))) %>%
  dplyr::select(id, SpeciesAuthor, MatrixPopulation, start,
                loglam, damp, gen, pmature, starts_with("id_"))


save(pt_shape_out, file = "analysis/full_pt_shape.RData")
save(pt_other_out, file = "analysis/full_pt_other.RData")

save(sd_shape_out, file = "analysis/full_sd_shape.RData")
save(sd_other_out, file = "analysis/full_sd_other.RData")



### testing new plot for appendix
dat_example <- pt_shape %>% 
  filter(Authors == "Martin; Meinke", MatrixPopulation == "Bull Flat")

stages <- c("Seedbank", "Seedling", "Small", "Medium", "Large")

N1 <- dat_example$rep_prop1[[1]]
matU <- dat_example$matU[[1]]
q <- dat_example$q[1]

df_lx <- tibble(
  lx = Rage::mpm_to_lx(matU, N1),
  hx = Rage::lx_to_hx(lx),
  x = seq_along(lx) - 1
) %>% filter(x <= 9)

proj <- popbio::pop.projection(matU, N1, iterations = 10)$stage.vectors %>% 
  t() %>% 
  apply(., 1, function(x) x / sum(x)) %>% 
  t() %>% 
  as.data.frame() %>% 
  setNames(stages) %>% 
  as_tibble() %>% 
  mutate(x = 1:n() - 1) %>% 
  tidyr::gather("stage", "tr", -x) %>% 
  mutate(stage = factor(stage, levels = stages))

p1 <- ggplot(proj, aes(x, tr, col = stage)) +
  geom_line(size = 2) +
  geom_vline(xintercept = q, linetype = 2) +
  scale_color_brewer(direction = 2, name = "Stage class") +
  scale_x_continuous(limits = c(0, 9), breaks = seq(0, 8, 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "Relative stage distribution") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank())

p2 <- ggplot(df_lx, aes(x, lx)) +
  geom_line() +
  geom_vline(xintercept = q, linetype = 2) +
  scale_x_continuous(limits = c(0, 9), breaks = seq(0, 8, 1)) +
  # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_log10(labels = formatC) +
  labs(x = NULL, y = "Survivorship") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank())

p3 <- ggplot(df_lx, aes(x, hx)) +
  geom_line() +
  geom_vline(xintercept = q, linetype = 2) +
  scale_x_continuous(limits = c(0, 9), breaks = seq(0, 8, 1)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
  # scale_y_log10(labels = formatC) +
  labs(x = "Time", y = "Mortality hazard") +
  theme_bw() +
  theme(panel.grid = element_blank())


library(patchwork)

g <- p1 / p2 / p3 + patchwork::plot_annotation(tag_levels = c("A", "B", "C"))

# print to screen
graphics.off()
quartz(height = 6, width = 6.5, dpi = 150); print(g)

# save png
# ggsave("img/appendix_qsd.png", g, height = 6, width = 6.5, units = "in", dpi = 600)

