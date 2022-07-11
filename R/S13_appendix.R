

### libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(patchwork)
library(popbio)
library(Rage)
library(tidyverse)
library(Rcompadre)
source("R/functions.R")


# rbeta(1, shape1 = 1 + k, shape2 = 1 + (n-k))
# rdirichlet(1, alpha = 1 + x)
# rgamma(1, shape = 1 + k, rate = n)


### Table S1-2
matU <- rbind(c(0.000, 0.000, 0.000, 0.000),
              c(0.857, 0.778, 0.000, 0.000),
              c(0.000, 0.056, 0.667, 0.277),
              c(0.000, 0.000, 0.222, 0.702))

matF <- rbind(c(0.000, 0.000, 0.000, 0.191),
              c(0.000, 0.000, 0.000, 0.000),
              c(0.000, 0.000, 0.000, 0.000),
              c(0.000, 0.000, 0.000, 0.000))

matA <- matU + matF

# stage-specific sample sizes
n <- c(14, 18, 9, 47)

# back-calculated transition frequencies
xU <- round(t(n * t(matU)))    # survival-related transitions
xF <- round(t(n * t(matF)))    # reproductive transitions
round(n * (1 - colSums(matU))) # deaths

# draw from sampling distribution
set.seed(987654321)
drawsU <- sim_U_wrapper(matU, N = n, nsim = 2000)
drawsF <- sim_U_wrapper(matF, N = n, nsim = 2000)
drawsA <- mapply(function(x, y) x + y, drawsU, drawsF, SIMPLIFY = F)

# MPM sampling distributions (draws 1, 2, ..., 2000)
lapply(drawsA[c(1, 2, 2000)], function(m) round(m, 3))

# derived param point estimates
round(lambda(matA), 2)
round(damping.ratio(matA), 2)
round(life_expect(matU), 2)

# derived param sampling distributions (draws 1, 2, ..., 2000)
round(sapply(drawsA[c(1, 2, 2000)], lambda), 2)
round(sapply(drawsA[c(1, 2, 2000)], damping.ratio), 2)
round(sapply(drawsU[c(1, 2, 2000)], life_expect), 2)



### Table S1-2

# component MPMs
mU1 <- rbind(c(0, 0), c(6/8, 5/14))
mU2 <- rbind(c(0, 0), c(5/8, 7/24))
mU3 <- rbind(c(0, 0), c(20/29, 3/23))

mF1 <- rbind(c(0, 14/14), c(0, 0))
mF2 <- rbind(c(0, 21/24), c(0, 0))
mF3 <- rbind(c(0, 24/23), c(0, 0))

mA1 <- mU1 + mF1
mA2 <- mU2 + mF2
mA3 <- mU3 + mF3

n1 <- c(8, 14)
n2 <- c(8, 24)
n3 <- c(29, 23)

# mean MPM
mU <- popbio::mean.list(list(mU1, mU2, mU3))
mA <- popbio::mean.list(list(mA1, mA2, mA3))

# point estimates for component and mean MPMs
round(mA, 2)
round(mA1, 2)
round(mA2, 2)
round(mA3, 2)

# point estimates for derived params from mean MPM
round(lambda(mA), 2)
round(damping.ratio(mA), 2)
round(life_expect(mU), 2)

# draw from the sampling distributions of component MPMs
set.seed(987654321)
drawsU1 <- sim_U_wrapper(mU1, N = n1, nsim = 2000)
drawsF1 <- sim_U_wrapper(mF1, N = n1, nsim = 2000)
drawsA1 <- mapply(function(x, y) x + y, drawsU, drawsF, SIMPLIFY = F)

drawsU2 <- sim_U_wrapper(mU2, N = n2, nsim = 2000)
drawsF2 <- sim_U_wrapper(mF2, N = n2, nsim = 2000)
drawsA2 <- mapply(function(x, y) x + y, drawsU, drawsF, SIMPLIFY = F)

drawsU3 <- sim_U_wrapper(mU3, N = n3, nsim = 2000)
drawsF3 <- sim_U_wrapper(mF3, N = n3, nsim = 2000)
drawsA3 <- mapply(function(x, y) x + y, drawsU, drawsF, SIMPLIFY = F)

# derive draws of mean MPM from draws of components
drawsU <- list(mean.list(list(drawsU1[[1]], drawsU2[[1]], drawsU3[[1]])),
               mean.list(list(drawsU1[[2]], drawsU2[[2]], drawsU3[[2]])),
               mean.list(list(drawsU1[[2000]], drawsU2[[2000]], drawsU3[[2000]])))

drawsF <- list(mean.list(list(drawsF1[[1]], drawsF2[[1]], drawsF3[[1]])),
               mean.list(list(drawsF1[[2]], drawsF2[[2]], drawsF3[[2]])),
               mean.list(list(drawsF1[[2000]], drawsF2[[2000]], drawsF3[[2000]])))

drawsA <- mapply(function(x, y) x + y, drawsU, drawsF, SIMPLIFY = F)


# draws from sampling distributions of component MPMs
lapply(drawsA1[c(1, 2, 2000)], function(x) round(x, 2))
lapply(drawsA2[c(1, 2, 2000)], function(x) round(x, 2))
lapply(drawsA3[c(1, 2, 2000)], function(x) round(x, 2))

# draws from sampling distribution of mean MPM
lapply(drawsA, function(m) round(m, 2))

# draws from sampling distribution of derived params from mean MPM
round(sapply(drawsA, lambda), 2)
round(sapply(drawsA, damping.ratio), 2)
round(sapply(drawsU, life_expect), 2)




### boundary estimates of survival
load("analysis/sd_scanga.RData")

scanga_t <- scanga_out %>% 
  filter(MatrixPopulation == "T")

(spp <- scanga_t$SpeciesAuthor[1])


matU <- scanga_t$mat[[1]]@matU
matU[,3][matU[,3] > 0] <- matU[,3][matU[,3] > 0] - (0.01 / 3)
matU[,4][matU[,4] > 0] <- matU[,4][matU[,4] > 0] - (0.01 / 4)
matU[,5][matU[,5] > 0] <- matU[,5][matU[,5] > 0] - (0.01 / 2)

9/14
round(matU, 2)

1 - colSums(matU) %>% round(2)

stages <- c("Seedling", "Small", "Medium", "Large", "X-Large")
stages <- factor(stages, levels = stages)

sigma_point <- scanga_t %>% 
  cdb_unnest() %>% 
  as_tibble() %>% 
  select(matU) %>% 
  mutate(n = 1:n(),
         sigma = map(matU, colSums),
         stage = list(stages)) %>% 
  select(-matU) %>% 
  unnest(c("sigma", "stage"))

sigma_sim <- scanga_t %>% 
  as_tibble() %>% 
  select(simU) %>% 
  unnest("simU") %>% 
  mutate(n = 1:n(),
         sigma = map(simU, colSums),
         stage = list(stages)) %>% 
  select(-simU) %>% 
  unnest(c("sigma", "stage"))


load("analysis/full_pt_shape.RData")
load("analysis/full_sd_shape.RData")

life_point <- pt_shape_out %>% 
  filter(SpeciesAuthor == "Trollius_laxus_2", MatrixPopulation == "T")

life_sim <- sd_shape_out %>% 
  filter(SpeciesAuthor == "Trollius_laxus_2", MatrixPopulation == "T")

mean(life_sim$L)
quantile(life_sim$L, c(0.5, 0.025, 0.975))


tt <- theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 11.5),
        axis.ticks = element_line(size = 0.4))

p1 <- ggplot(sigma_sim, aes(sigma)) +
  geom_density(fill = "darkred", color = NA, alpha = 0.5) +
  geom_vline(data = sigma_point, aes(xintercept = sigma), linetype = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 1)) +
  coord_cartesian(ylim = c(0, 20)) +
  facet_wrap(~ stage, nrow = 1) +
  labs(x = "Stage-specific survival probability",
       y = "Density") +
  tt

p2 <- ggplot(life_sim, aes(L)) +
  geom_density(fill = "darkred", color = NA, alpha = 0.5) +
  geom_vline(data = life_point, aes(xintercept = L_pt), linetype = 2) +
  scale_x_continuous(limits = c(0, 148)) +
  labs(x = "Mature life expectancy (years)",
       y = "Density") +
  tt

g1 <- p1 / p2 + plot_layout(heights = c(0.7, 1))


dev.off()
quartz(height = 4, width = 6.25, dpi = 160)
print(g1)

ggsave("img/boundary.png", g1, height = 4, width = 6.25, units = "in", dpi = 300)



scanga_traj_pt <- pt_shape %>% 
  filter(SpeciesAuthor == spp, MatrixPopulation == "T") %>% 
  mutate(rep = 1:n()) %>% 
  select(rep, lx) %>% 
  unnest(cols = "lx") %>% 
  group_by(rep) %>% 
  mutate(x = seq_along(lx),
         hx = Rage::lx_to_hx(lx)) %>% 
  ungroup()

scanga_traj_sim <- sd_shape %>% 
  filter(SpeciesAuthor == spp, MatrixPopulation == "T") %>% 
  mutate(rep = 1:n()) %>% 
  select(rep, lx) %>% 
  unnest(cols = "lx") %>% 
  group_by(rep) %>% 
  mutate(x = seq_along(lx),
         hx = Rage::lx_to_hx(lx)) %>% 
  ungroup() # %>% 
  # group_by(x) %>% 
  # summarize(lx_low = quantile(lx, 0.025),
  #           lx_upp = quantile(lx, 0.975))

p3 <- ggplot(scanga_traj_sim, aes(x, lx)) +
  geom_line(aes(group = rep), color = "darkred", alpha = 0.025) +
  # geom_ribbon(aes(ymin = lx_low, ymax = lx_upp), fill = "darkred", alpha = 0.5) +
  geom_line(data = scanga_traj_pt, linetype = 2) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  # scale_y_log10() +
  labs(x = "Age from reproductive maturity (years)",
       y = expression(paste("Survivorship (", italic(l[x]), ")"))) +
  tt


g2 <- p1 / p2 / p3 +
  plot_layout(heights = c(0.8, 1, 1)) +
  plot_annotation(tag_levels = "a")


dev.off()
quartz(height = 6, width = 6.25, dpi = 160)
print(g2)

ggsave("img/boundary2.png", g2, height = 6, width = 6.25, units = "in", dpi = 300)



### illustrate pooled

### Kiviniemi
spp <- "Agrimonia_eupatoria"

compadre %>% 
  filter(SpeciesAuthor == spp, MatrixPopulation == "A") %>% 
  filter(MatrixTreatment == 'Unmanipulated') %>%
  cdb_glimpse("MatrixComposite")

kiviniemi_n <- read_csv("data/studies/kiviniemi_n.csv") %>% 
  filter(SpeciesAccepted == "Agrimonia eupatoria", MatrixPopulation == "A") %>% 
  group_by(SpeciesAccepted, MatrixPopulation, MatrixStartYear) %>% 
  summarize(N = list(N)) %>% 
  ungroup()

kiviniemi <- compadre %>% 
  filter(SpeciesAuthor == spp, MatrixPopulation == "A") %>% 
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
  filter(SpeciesAuthor == spp, MatrixPopulation == "A") %>% 
  filter(MatrixComposite == "Mean") %>% 
  filter(MatrixTreatment == "Unmanipulated") %>% 
  left_join(npool)

