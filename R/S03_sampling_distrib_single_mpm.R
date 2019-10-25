
### libraries
library(tidyverse)
library(Rcompadre)
library(Rage)
library(popbio)
library(gridExtra)
source("R/functions.R")


### load compadre
compadre <- cdb_fetch("data/COMPADRE_v.X.X.X_Corrected.RData")


### kiviniemi example
kiviniemi_n <- read_csv("data/studies/kiviniemi_n.csv") %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  summarize(N = list(N)) %>% 
  ungroup()

# focal species and population
kiviniemi <- compadre %>%
  filter(SpeciesAuthor == "Agrimonia_eupatoria",
         MatrixPopulation == "B",
         MatrixComposite == "Individual",
         MatrixStartYear == 1997) %>% 
  cdb_unnest() %>% 
  left_join(kiviniemi_n, by = c("MatrixPopulation", "MatrixStartYear"))

# short form stage class names
kiviniemi$MatrixClassAuthor[[1]]
stage_names <- c('Sdl.', 'Juv.', 'Veg.', 'Rep.')

# stage sample-sizes for focal year and population
kiviniemi_n1 <- kiviniemi_n %>% 
  filter(MatrixPopulation %in% kiviniemi$MatrixPopulation,
         MatrixStartYear %in% kiviniemi$MatrixStartYear) %>% 
  unnest(cols = "N") %>% 
  mutate(from_col = 1:n())

# convert mpm to flat form
df_mpm <- mpm_flatten(kiviniemi$matA[[1]],
                      kiviniemi$matU[[1]],
                      kiviniemi$matF[[1]],
                      kiviniemi$matC[[1]],
                      stage_names)

# arrange data for plotting
df_plot <- df_mpm %>% 
  left_join(kiviniemi_n1, by = "from_col") %>% 
  mutate(num = round(A * N, 1)) %>% 
  mutate(n_lab = ifelse(A == 0, "", paste0(num, "/", N))) %>%
  mutate(A_num = ifelse(A == 0, NA_real_, A)) %>% 
  mutate(A = ifelse(A == 0, "", sprintf("%.2f", A))) %>% 
  mutate(U = ifelse(U == 0, "", sprintf("%.2f", U))) %>% 
  mutate(F = ifelse(F == 0, "", sprintf("%.2f", F))) %>% 
  mutate(C = ifelse(C == 0, "", sprintf("%.2f", C)))

df_sdist <- df_mpm %>% 
  left_join(kiviniemi_n1, by = "from_col") %>% 
  mutate(x = round(A * N, 1)) %>% 
  mutate(x = ifelse(A == 0, NA, x)) %>% 
  mutate(fec = ifelse(F == 0, FALSE, TRUE)) %>% 
  group_by(to_col, from_col) %>% 
  do(dens_fn(.$x, .$N, .$fec)) %>% 
  ungroup() %>% 
  left_join(df_mpm, by = c("to_col", "from_col"))

df_pt <- df_sdist %>% 
  group_by(to_name, from_name) %>% 
  filter(!is.na(pp), pp == max(pp)) %>% 
  ungroup()

df_rect <- df_mpm %>% 
  mutate(x1 = ifelse(A == 0, -Inf, NA),
         x2 = ifelse(A == 0, Inf, NA)) %>% 
  mutate(y1 = x1, y2 = x2)


### Figure 1 (top): MPM components
p1a <- ggplot(df_plot) +
  geom_text(aes(label = A, x = 1, y = 1), size = 3.2) +
  facet_grid(to_name ~ from_name, switch = "y") +
  labs(x = NULL, y = NULL) +
  ggtitle("In COMPADRE") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0),
        text = element_text(size = 11.7),
        plot.margin = margin(5, 16, 2, 2),
        panel.spacing = unit(1.3, "pt"))

p1b <- ggplot(df_plot) +
  geom_text(aes(label = n_lab, x = 1, y = 1), size = 3.2) +
  facet_grid(to_name ~ from_name, switch = "y") +
  labs(x = NULL, y = NULL) +
  ggtitle("Raw data") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0),
        text = element_text(size = 11.7),
        plot.margin = margin(5, 16, 2, 2),
        panel.spacing = unit(1.3, "pt"))

p1c <- ggplot(df_sdist) +
  geom_ribbon(aes(x = p, ymin = 0, ymax = pp), fill = "darkred", alpha = 0.5) +
  geom_segment(data = df_pt, aes(x = p, y = 0, xend = p, yend = pp + 0.1), size = 0.3) +
  geom_rect(data = df_rect, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            fill = "grey93") +
  facet_grid(to_name ~ from_name, switch = "y") +
  coord_flip() +
  scale_y_reverse() +
  scale_x_continuous(limits = c(0, 1),
                     expand = c(0.06, 0),
                     position = "top",
                     breaks = seq(0, 1, 0.2),
                     labels = c(" ", " ", " ", " ", " ", " ")) +
  labs(x = NULL, y = NULL) +
  ggtitle("Sampling distributions") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        # axis.ticks.y = element_line(size = 0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.15),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0),
        text = element_text(size = 11.7),
        plot.margin = margin(5, 2, 2, 2),
        panel.spacing = unit(1.3, "pt"))

p1 <- cbind(ggplotGrob(p1a), ggplotGrob(p1b), ggplotGrob(p1c), size = "first")

dev.off()
quartz(height = 2.2, width = 6.5, dpi = 140)
grid.arrange(p1)

# ggsave("img/raw/Fig_1a.png", p1, height = 2.2, width = 6.5, units = "in", dpi = 300)



### Figure 1 (bottom): Derived parameters

# possible transitions
posU <- mean(kiviniemi$matU) > 0
posF <- mean(kiviniemi$matF) > 0
posA <- mean(kiviniemi$matA) > 0

# generate MPMs based on sampling ditributions of transitions rates
individ_sim <- kiviniemi %>%
  mutate(simU = map2(matU, N, ~ sim_U_wrapper(.x, posU, .y, 2000))) %>%
  mutate(simF = map2(matF, N, ~ sim_F_wrapper(.x, posF, .y, 2000))) %>%
  mutate(simA = map2(simF, simU, ~ map2(.x, .y, `+`))) %>%
  as_tibble() %>% 
  select(simU, simF, simA) %>%
  unnest() %>%
  mutate(rep = 1:n()) %>%
  select(rep, simA, simU, simF)

# paramater factor levels
par_lev <- c("lambda", "rho", "v[3]", "E[list(4,4)]", "T", "l[0]")
par_lab <- paste0("italic(", par_lev, ")")

# derived parameters, sampling distributions
deriv_param <- individ_sim %>% 
  mutate(lambda = map_dbl(simA, lambda)) %>% 
  mutate(rho = map_dbl(simA, damping.ratio)) %>% 
  mutate(`v[3]` = map_dbl(simA, ~ reproductive.value(.x)[3])) %>%
  mutate(`E[list(4,4)]` = map_dbl(simA, ~ elasticity(.x)[4,4])) %>% 
  mutate(`R[0]` = map2_dbl(simU, simF, net_repro_rate)) %>%
  mutate(`T` = map2_dbl(simU, simF, gen_time)) %>% 
  mutate(`l[0]` = map_dbl(simU, life_expect)) %>%
  select(-matches("sim"), -`R[0]`) %>% 
  gather(par, value, lambda:`l[0]`) %>% 
  mutate(par = factor(par, levels = par_lev, labels = par_lab))

# derived parameters, point estimates
deriv_pt <- kiviniemi %>% 
  as_tibble() %>% 
  select(matA, matU, matF) %>% 
  mutate(lambda = map_dbl(matA, lambda)) %>% 
  mutate(rho = map_dbl(matA, damping.ratio)) %>% 
  mutate(`v[3]` = map_dbl(matA, ~ reproductive.value(.x)[3])) %>%
  mutate(`E[list(4,4)]` = map_dbl(matA, ~ elasticity(.x)[4,4])) %>% 
  mutate(`R[0]` = map2_dbl(matU, matF, net_repro_rate)) %>%
  mutate(`T` = map2_dbl(matU, matF, gen_time)) %>% 
  mutate(`l[0]` = map_dbl(matU, life_expect)) %>%
  select(-starts_with("mat"), -`R[0]`) %>% 
  gather(par, value, lambda:`l[0]`) %>% 
  mutate(par = factor(par, levels = par_lev, labels = par_lab))

# plot
deriv_plot <- deriv_param %>% 
  filter(!(par == "italic(v[3])" & value > 18)) %>% 
  filter(!(par == "italic(T)" & value > 90))

p2 <- ggplot(deriv_plot) +
  geom_density(aes(value), fill = "darkred", alpha = 0.5, size = 0) +
  geom_vline(data = deriv_pt, aes(xintercept = value)) +
  facet_wrap(~ par, scales = "free", labeller = label_parsed, nrow = 1) +
  labs(x = "Parameter estimate", y = "Probability density") +
  ggtitle("Derived parameters") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.3),
        plot.title = element_text(hjust = 0, face = "bold", vjust = 0, size = 12),
        text = element_text(size = 11.7),
        strip.text = element_text(size = 12))


dev.off()
quartz(height = 2, width = 6.5, dpi = 120)
print(p2)

# ggsave("img/raw/Fig_1b.png", p2, height = 2, width = 6.5, units = "in", dpi = 300)



### table of posterior quantiles for derived parameters
deriv_param %>% 
  group_by(par) %>% 
  summarize(med = quantile(value, 0.500),
            low = quantile(value, 0.025),
            upp = quantile(value, 0.975)) %>% 
  left_join(deriv_pt, by = "par")

