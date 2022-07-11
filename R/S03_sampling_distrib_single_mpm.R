
### libraries
library(tidyverse)
library(Rcompadre)
library(Rage)
library(popbio)
library(gridExtra)
source("R/functions.R")


### load compadre
compadre <- Rcompadre::cdb_fetch("data/COMPADRE_v.X.X.X_Corrected.RData")


### kiviniemi example
kiviniemi_n <- readr::read_csv("data/studies/kiviniemi_n.csv", show_col_types = FALSE) %>% 
  group_by(MatrixPopulation, MatrixStartYear) %>% 
  summarize(N = list(N), .groups = "drop")

# focal species and population
kiviniemi <- compadre %>%
  cdb_unnest() %>% 
  as_tibble() %>% 
  filter(SpeciesAuthor == "Agrimonia_eupatoria",
         MatrixPopulation == "B",
         MatrixComposite == "Individual",
         MatrixStartYear == 1997) %>% 
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
  theme_bw() +
  theme(
    strip.background = element_rect(color = "grey80", fill = "grey85", size = 0.4),
    strip.text = element_text(margin = margin(0.2, 0.25, 0.2, 0.25, "lines")),
    # strip.background = element_blank(),
    panel.border = element_rect(color = "grey80"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11.5, vjust = 0),
    text = element_text(size = 11.7),
    plot.margin = margin(5, 16, 2, 2),
    panel.spacing = unit(1.3, "pt")
  )

p1b <- ggplot(df_plot) +
  geom_text(aes(label = n_lab, x = 1, y = 1), size = 3.2) +
  facet_grid(to_name ~ from_name, switch = "y") +
  labs(x = NULL, y = NULL) +
  ggtitle("Raw data") +
  theme_bw() +
  theme(
    strip.background = element_rect(color = "grey80", fill = "grey85", size = 0.4),
    strip.text = element_text(margin = margin(0.2, 0.25, 0.2, 0.25, "lines")),
    panel.border = element_rect(color = "grey80"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11.5, vjust = 0),
    text = element_text(size = 11.7),
    plot.margin = margin(5, 16, 2, 2),
    panel.spacing = unit(1.3, "pt")
  )

p1c <- ggplot(df_sdist) +
  geom_ribbon(aes(x = p, ymin = 0, ymax = pp), fill = "darkred", alpha = 0.4) +
  geom_segment(data = df_pt, aes(x = p, y = 0, xend = p, yend = pp + 0.1), size = 0.3) +
  geom_rect(data = df_rect, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            fill = NA) +
  facet_grid(to_name ~ from_name, switch = "y") +
  coord_flip() +
  scale_y_reverse() +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 1),
    expand = c(0, 0.11),
    position = "top"
  ) +
  labs(x = NULL, y = NULL) +
  ggtitle("Sampling distributions") +
  theme_bw() +
  theme(
    strip.background = element_rect(color = "grey80", fill = "grey85", size = 0.4),
    strip.text = element_text(margin = margin(0.2, 0.25, 0.2, 0.25, "lines")),
    panel.border = element_rect(color = "grey80"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 0.2, color = "grey80"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11.5, vjust = 0),
    text = element_text(size = 11.7),
    plot.margin = margin(5, 2, 2, 2),
    panel.spacing = unit(1.3, "pt")
  )

p1 <- cbind(ggplotGrob(p1a), ggplotGrob(p1b), ggplotGrob(p1c), size = "first")

graphics.off()
quartz(height = 2.2, width = 6.5, dpi = 240)
grid.arrange(p1)

# ggsave("img/raw/Fig_1a.png", p1, height = 2.2, width = 6.5, units = "in", dpi = 400)



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
  unnest(cols = c(simU, simF, simA)) %>%
  mutate(rep = 1:n()) %>%
  select(rep, simA, simU, simF)


# paramater factor levels
par_lev <- c("lambda", "rho", "v[3]", "E[list(4,4)]", "T", "l[0]")
par_lab <- paste0("italic(", par_lev, ")")

par_lab <- c(
  "Pop.~growth~rate~(italic(lambda))",
  "Damping~ratio~(italic(rho))",
  "Repro.~value~(Veg.)",
  "Elasticity",
  "Generation~time~(italic(T))",
  "Life~expectancy~(italic(l[0]))"
)


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
  mutate(par = factor(par, levels = par_lev, labels = par_lab)) %>% 
  filter(!par %in% c("Repro.~value~(Veg.)", "Elasticity"))

# plot
deriv_plot <- deriv_param %>% 
  filter(!par %in% c("Repro.~value~(Veg.)", "Elasticity")) %>% 
  filter(!(par == "Generation~time~(italic(T))" & value > 90)) %>%
  filter(!(par == "Life~expectancy~(italic(l[0]))" & value > 35))

p2 <- ggplot(deriv_plot) +
  geom_density(aes(value), fill = "darkred", alpha = 0.4, size = 0) +
  geom_vline(data = deriv_pt, aes(xintercept = value)) +
  facet_wrap(~ par, scales = "free", labeller = label_parsed, nrow = 1) +
  labs(x = "Parameter estimate", y = "Prob. density") +
  ggtitle("Derived parameters") +
  theme_bw() +
  theme(
    strip.background = element_rect(color = "grey80", fill = "grey85", size = 0.4),
    panel.border = element_rect(color = "grey80"),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_line(size = 0.4, color = "grey80"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(size = 0.3, color = "grey80"),
    plot.title = element_text(hjust = 0, face = "bold", vjust = 0, size = 12),
    text = element_text(size = 11.7),
    strip.text = element_text(size = 9, margin = margin(0.15, 0, 0.15, 0, "lines"))
  )


graphics.off()
quartz(height = 2, width = 6.5, dpi = 220)
print(p2)

# ggsave("img/raw/Fig_1b.png", p2, height = 2, width = 6.5, units = "in", dpi = 300)



### test approach based on patchwork
# library(patchwork)
# 
# plotfn <- function(var, title = str2lang(var)) {
#   deriv_param_foc <- deriv_param %>%
#     filter(par == var)
# 
#   deriv_pt_foc <- deriv_pt %>%
#     filter(par == var)
#   
#   ggplot(deriv_param_foc) +
#     geom_density(aes(value), fill = "darkred", alpha = 0.4, size = 0) +
#     geom_vline(data = deriv_pt_foc, aes(xintercept = value)) +
#     # labs(x = "Parameter estimate", y = "Probability density") +
#     labs(x = NULL, y = NULL) +
#     # ggtitle(title) +
#     theme_bw() +
#     theme(
#       strip.background = element_rect(color = "grey80", fill = "grey85", size = 0.4),
#       panel.border = element_rect(color = "grey80"),
#       panel.grid = element_blank(),
#       axis.text.y = element_blank(),
#       axis.ticks = element_line(size = 0.4),
#       axis.ticks.y = element_blank(),
#       axis.ticks.x = element_line(size = 0.3),
#       text = element_text(size = 11),
#       plot.title = element_text(hjust = 0, face = "bold", vjust = 0, size = 10)
#   )
# }
# 
# p2a <- plotfn("italic(lambda)") +
#   scale_x_continuous(limits = c(0.74, 1.1)) +
#   # labs(x = NULL, y = "Probability density") +
#   ggtitle(expression(paste("Popn. growth rate (", italic(lambda), ")")))
# 
# p2b <- plotfn("italic(rho)") +
#   scale_x_continuous(limits = c(0.8, 2.8)) +
#   ggtitle(expression(paste("Damping ratio (", italic(rho), ")")))
# 
# p2c <- plotfn("italic(v[3])") +
#   scale_x_log10(
#     limits = c(0.1, 100),
#     breaks = c(0.1, 1, 10, 100),
#     labels = as.character(c(0.1, 1, 10, 100))
#   ) +
#   labs(x = "Parameter estimate", y = NULL) +
#   coord_cartesian(clip = "off") +
#   ggtitle(expression(paste("Repro value (Veg.)")))
# 
# p2d <- plotfn("italic(T)") +
#   scale_x_continuous(limits = c(0, 100), breaks = seq(0, 90, 30)) +
#   ggtitle(expression(paste("Generation time (", italic(T), ")")))
# 
# p2e <- plotfn("italic(l[0])") +
#   scale_x_continuous(limits = c(0, 30)) +
#   ggtitle(expression(paste("Life expectancy (", italic(l[0]), ")")))
# 
# p1 <- p1a + p1b + p1c + plot_layout(nrow = 1)
# p2 <- p2a + p2b + p2c + p2d + p2e + plot_layout(nrow = 1)
# 
# g <- (p1a + p1b + p1c) / p2 +
#   plot_layout(heights = c(0.6, 0.4))
# 
# g <- p1 / wrap_elements(full = p2) + plot_layout(heights = c(0.65, 0.35))
# 
# 
# graphics.off()
# quartz(height = 4, width = 6.5, dpi = 220)
# print(g)
# 
# top_row <- plot_grid(
#   p1a, p1b, p1c, labels = c('A', 'B', 'C'),
#   nrow = 1, label_size = 12
# )
# 
# g <- plot_grid(
#   top_row, p2,
#   labels = c("", "D"),
#   label_size = 12, ncol = 1,
#   rel_heights = c(1.2, 1)
# )
# 
# graphics.off()
# quartz(height = 4, width = 6.5, dpi = 220)
# print(g)
# 
# ggsave("img/fig_1_new.png", g, height = 4, width = 6.5, units = "in", dpi = 600)



### table of posterior quantiles for derived parameters
deriv_param %>% 
  group_by(par) %>% 
  summarize(med = quantile(value, 0.500),
            low = quantile(value, 0.025),
            upp = quantile(value, 0.975)) %>% 
  left_join(deriv_pt, by = "par")

