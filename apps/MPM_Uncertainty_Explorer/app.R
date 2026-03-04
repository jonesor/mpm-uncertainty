# MPM Uncertainty Explorer (tabbed dimensions: 2x2, 3x3, 4x4)

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

set.seed(5654)

quantity_levels <- c(
  "Lambda",
  "Mature life expectancy (L)",
  "Generation time (T)",
  "Damping ratio (rho)"
)

u_col <- "#35b779"
f_col <- "#e78ac3"
point_col <- "#c7352a"
post_col <- "#7a7a7a"
derived_col <- "#e68613"

fmt_val <- function(x) formatC(x, format = "f", digits = 2, drop0trailing = TRUE, decimal.mark = ".")

snap_u_count_column <- function(k_col, n_col) {
  k <- suppressWarnings(as.integer(round(k_col)))
  k[!is.finite(k)] <- 0L
  k <- pmax(k, 0L)
  overflow <- sum(k) - n_col
  if (overflow > 0) {
    ord <- order(k, decreasing = TRUE)
    idx_ptr <- 1L
    max_iter <- 4000L
    iter <- 0L
    while (overflow > 0 && iter < max_iter) {
      idx <- ord[((idx_ptr - 1L) %% length(ord)) + 1L]
      if (!is.na(idx) && k[idx] > 0L) {
        k[idx] <- k[idx] - 1L
        overflow <- overflow - 1L
      }
      idx_ptr <- idx_ptr + 1L
      iter <- iter + 1L
    }
  }
  k
}

build_defaults <- function(k) {
  if (k == 3) {
    n <- c(20L, 15L, 10L)
    U <- matrix(0L, 3, 3)
    U[1, 1] <- 5L
    U[2, 1] <- 7L
    U[2, 2] <- 8L
    U[3, 2] <- 5L
    U[3, 3] <- 7L
    F <- matrix(0L, 3, 3)
    F[1, 3] <- 20L
    return(list(n = n, U = U, F = F))
  }
  if (k == 2) {
    n <- c(20L, 15L)
    U <- matrix(c(8L, 0L,
                  6L, 7L), nrow = 2, byrow = TRUE)
    F <- matrix(0L, 2, 2)
    F[1, 2] <- 15L
    return(list(n = n, U = U, F = F))
  }
  # k == 4
  n <- c(20L, 15L, 10L, 8L)
  U <- matrix(0L, 4, 4)
  U[1, 1] <- 6L
  U[2, 1] <- 7L
  U[2, 2] <- 6L
  U[3, 2] <- 5L
  U[3, 3] <- 4L
  U[4, 3] <- 3L
  U[4, 4] <- 5L
  F <- matrix(0L, 4, 4)
  F[1, 4] <- 16L
  list(n = n, U = U, F = F)
}

calc_density_beta <- function(p_hat, n) {
  if (length(n) != 1 || is.na(n) || !is.finite(n) || n <= 0) {
    return(tibble(x = seq(0, 1, by = 0.01), d = NA_real_))
  }
  if (length(p_hat) != 1 || is.na(p_hat) || !is.finite(p_hat) || p_hat < 0 || p_hat > 1) {
    return(tibble(x = seq(0, 1, by = 0.01), d = NA_real_))
  }
  x_obs <- round(p_hat * n)
  grid <- seq(0, 1, by = 0.01)
  dens <- dbeta(grid, x_obs + 1, n - x_obs + 1)
  tibble(x = grid, d = dens / max(dens, na.rm = TRUE))
}

calc_density_gamma <- function(rate_hat, n) {
  if (length(n) != 1 || is.na(n) || !is.finite(n) || n <= 0) {
    return(tibble(x = seq(0, 12, length.out = 400), d = NA_real_))
  }
  if (length(rate_hat) != 1 || is.na(rate_hat) || !is.finite(rate_hat) || rate_hat < 0) {
    return(tibble(x = seq(0, 12, length.out = 400), d = NA_real_))
  }
  y_obs <- round(rate_hat * n)
  x_max <- max(12, qgamma(0.995, shape = y_obs + 1, rate = n))
  grid <- seq(0, x_max, length.out = 400)
  dens <- dgamma(grid, shape = y_obs + 1, rate = n)
  tibble(x = grid, d = dens / max(dens, na.rm = TRUE))
}

rdirichlet1 <- function(alpha) {
  z <- rgamma(length(alpha), shape = alpha, rate = 1)
  z / sum(z)
}

spectral_summary <- function(A) {
  eig <- tryCatch(eigen(A)$values, error = function(e) rep(NA_complex_, nrow(A)))
  mod_vals <- sort(Mod(eig), decreasing = TRUE)
  lambda1 <- if (all(is.na(mod_vals))) NA_real_ else Re(eig[which.max(Mod(eig))])
  rho <- if (length(mod_vals) >= 2 && is.finite(mod_vals[2]) && mod_vals[2] > 0) mod_vals[1] / mod_vals[2] else NA_real_
  list(lambda = lambda1, damping = rho)
}

life_expectancy_from_start <- function(U, start_stage) {
  if (any(is.na(U)) || is.na(start_stage) || start_stage < 1 || start_stage > nrow(U)) return(NA_real_)
  eig <- tryCatch(eigen(U)$values, error = function(e) NA_complex_)
  rad <- suppressWarnings(max(Mod(eig), na.rm = TRUE))
  if (!is.finite(rad) || rad >= 1) return(NA_real_)
  inv <- tryCatch(solve(diag(nrow(U)) - U), error = function(e) NULL)
  if (is.null(inv)) return(NA_real_)
  as.numeric(colSums(inv)[start_stage])
}

generation_time_safe <- function(U, F) {
  if (!requireNamespace("Rage", quietly = TRUE)) return(NA_real_)
  out <- tryCatch(Rage::gen_time(U, F), error = function(e) NA_real_)
  as.numeric(out)
}

theme_app_plot <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#dce6e0", linewidth = 0.3),
      strip.background = element_rect(fill = "#eef4f0", color = "#c9d6cf", linewidth = 0.5),
      strip.text = element_text(face = "bold", color = "#1f2a24", size = 8.5),
      plot.title = element_text(face = "bold", hjust = 0.5, color = "#1f2a24"),
      plot.subtitle = element_text(hjust = 0.5, color = "#4d6358", size = 9.5),
      axis.title = element_text(color = "#24342d"),
      axis.text = element_text(color = "#2c3f36")
    )
}

id_n <- function(prefix, j) paste0(prefix, "_n_", j)
id_u <- function(prefix, i, j) paste0(prefix, "_u_", i, "_", j)
id_f <- function(prefix, i, j) paste0(prefix, "_f_", i, "_", j)

panel_inputs <- function(prefix, k, defs) {
  tagList(
    tags$h4("Stage sample sizes"),
    fluidRow(lapply(seq_len(k), function(j) {
      column(width = max(2, floor(12 / k)), numericInput(id_n(prefix, j), paste0("n[", j, "]"), value = defs$n[j], min = 1, step = 1))
    })),
    tags$h4("U matrix"),
    tags$small("Counts of surviving/transitioning individuals."),
    fluidRow(lapply(seq_len(k), function(j) column(width = max(2, floor(12 / k)), tags$b(paste("from", j))))),
    lapply(seq_len(k), function(i) {
      fluidRow(
        column(12, tags$small(tags$b(paste("to", i)))),
        lapply(seq_len(k), function(j) {
          column(width = max(2, floor(12 / k)), numericInput(id_u(prefix, i, j), paste0("u[", i, ",", j, "]"), value = defs$U[i, j], min = 0, step = 1))
        })
      )
    }),
    tags$h4("F matrix"),
    tags$small("Counts of recruits from each source stage."),
    tags$p("Note: each F entry is a total recruit count across all individuals in the source stage.", style = "color: #c00000; font-weight: 600;"),
    fluidRow(lapply(seq_len(k), function(j) column(width = max(2, floor(12 / k)), tags$b(paste("from", j))))),
    lapply(seq_len(k), function(i) {
      fluidRow(
        column(12, tags$small(tags$b(paste("to", i)))),
        lapply(seq_len(k), function(j) {
          column(width = max(2, floor(12 / k)), numericInput(id_f(prefix, i, j), paste0("f[", i, ",", j, "]"), value = defs$F[i, j], min = 0, step = 1))
        })
      )
    })
  )
}

ui <- fluidPage(
  tags$head(tags$style(HTML(" 
      .well { border-radius: 14px; border: 1px solid #d7e3dd; background: linear-gradient(180deg, #ffffff 0%, #f8fcfa 100%); box-shadow: 0 8px 24px rgba(31, 42, 36, 0.08); }
      .plot-card { border: 1px solid #d7e3dd; border-radius: 14px; background: #ffffff; padding: 8px 10px 2px 10px; margin-bottom: 10px; box-shadow: 0 8px 24px rgba(31, 42, 36, 0.06); }
      .table-card { border: 1px solid #d7e3dd; border-radius: 14px; background: #ffffff; padding: 10px 14px; box-shadow: 0 8px 24px rgba(31, 42, 36, 0.06); }
      .details-overlay { position: fixed; inset: 0; background: rgba(15, 25, 20, 0.22); z-index: 2500; }
      .details-panel { position: fixed; top: 0; right: 0; width: 430px; max-width: 94vw; height: 100vh; overflow-y: auto; background: #ffffff; border-left: 1px solid #d7e3dd; box-shadow: -10px 0 30px rgba(31, 42, 36, 0.18); z-index: 2600; padding: 14px 16px 18px 16px; }
      .details-fab-wrap { position: fixed; right: 18px; bottom: 18px; z-index: 2400; }
      .details-fab-wrap .btn { border-radius: 999px; padding: 10px 16px; font-weight: 600; box-shadow: 0 8px 20px rgba(31, 42, 36, 0.2); }
    "))),
  titlePanel("MPM Uncertainty Explorer"),
  sidebarLayout(
    sidebarPanel(
      tags$p("Choose a matrix size, enter stage sample sizes and observed counts for U and F, then compare point estimates with posterior uncertainty in the plots and summary table."),
      tabsetPanel(
        id = "dim_tab",
        tabPanel("2x2", panel_inputs("d2", 2, build_defaults(2))),
        tabPanel("3x3", panel_inputs("d3", 3, build_defaults(3))),
        tabPanel("4x4", panel_inputs("d4", 4, build_defaults(4)))
      ),
      checkboxInput("fix_zero_structural", "Treat entered zeros as structural zeros (fixed at 0)", value = TRUE),
      selectInput("nsim", "Posterior draws", choices = c("300", "500", "1000"), selected = "500"),
      actionButton("reset_defaults", "Reset active tab defaults")
    ),
    mainPanel(
      fluidRow(
        column(8,
               div(class = "plot-card", plotOutput("u_plot", height = "380px")),
               div(class = "plot-card", plotOutput("f_plot", height = "380px"))
        ),
        column(4, div(class = "plot-card", plotOutput("derived_plot", height = "800px")))
      ),
      tags$hr(),
      div(class = "table-card",
          tags$h4("Derived Quantity Summary", style = "margin-top: 0; margin-bottom: 10px;"),
          tableOutput("summary_table")
      )
    )
  ),
  div(class = "details-fab-wrap", actionButton("toggle_details", "Details")),
  uiOutput("details_ui")
)

server <- function(input, output, session) {
  details_open <- reactiveVal(FALSE)
  observeEvent(input$toggle_details, details_open(TRUE))
  observeEvent(input$close_details, details_open(FALSE))

  active_cfg <- reactive({
    switch(input$dim_tab,
           "2x2" = list(prefix = "d2", k = 2L, defs = build_defaults(2)),
           "3x3" = list(prefix = "d3", k = 3L, defs = build_defaults(3)),
           "4x4" = list(prefix = "d4", k = 4L, defs = build_defaults(4)),
           list(prefix = "d3", k = 3L, defs = build_defaults(3)))
  })

  observeEvent(input$reset_defaults, {
    cfg <- active_cfg()
    for (j in seq_len(cfg$k)) updateNumericInput(session, id_n(cfg$prefix, j), value = cfg$defs$n[j])
    for (i in seq_len(cfg$k)) {
      for (j in seq_len(cfg$k)) {
        updateNumericInput(session, id_u(cfg$prefix, i, j), value = cfg$defs$U[i, j])
        updateNumericInput(session, id_f(cfg$prefix, i, j), value = cfg$defs$F[i, j])
      }
    }
    updateSelectInput(session, "nsim", selected = "500")
    updateCheckboxInput(session, "fix_zero_structural", value = TRUE)
  })

  inputs <- reactive({
    cfg <- active_cfg()
    n_vec <- map_dbl(seq_len(cfg$k), function(j) {
      v <- input[[id_n(cfg$prefix, j)]]
      if (is.null(v) || !is.finite(v)) cfg$defs$n[j] else as.numeric(v)
    }) %>% round() %>% pmax(1)

    U <- matrix(0, nrow = cfg$k, ncol = cfg$k)
    F <- matrix(0, nrow = cfg$k, ncol = cfg$k)
    for (i in seq_len(cfg$k)) {
      for (j in seq_len(cfg$k)) {
        uij <- suppressWarnings(as.numeric(input[[id_u(cfg$prefix, i, j)]]))
        fij <- suppressWarnings(as.numeric(input[[id_f(cfg$prefix, i, j)]]))
        U[i, j] <- if (is.null(uij) || !is.finite(uij)) cfg$defs$U[i, j] else uij
        F[i, j] <- if (is.null(fij) || !is.finite(fij)) cfg$defs$F[i, j] else fij
      }
    }
    list(k = cfg$k, n = n_vec, U = U, F = F)
  })

  validity <- reactive({
    x <- inputs()
    u_ok <- all(is.finite(x$U)) && all(x$U >= 0)
    f_ok <- all(is.finite(x$F)) && all(x$F >= 0)
    n_ok <- all(is.finite(x$n)) && all(x$n > 0)
    list(ok = u_ok && f_ok && n_ok)
  })

  observed_data <- reactive({
    req(validity()$ok)
    x <- inputs()
    U_counts <- matrix(0, nrow = x$k, ncol = x$k)
    F_counts <- matrix(0, nrow = x$k, ncol = x$k)
    for (j in seq_len(x$k)) U_counts[, j] <- snap_u_count_column(x$U[, j], x$n[j])
    for (j in seq_len(x$k)) for (i in seq_len(x$k)) F_counts[i, j] <- pmax(as.integer(round(x$F[i, j])), 0L)
    U_hat <- sweep(U_counts, 2, x$n, "/")
    F_hat <- sweep(F_counts, 2, x$n, "/")
    list(k = x$k, n = x$n, U_true = x$U, F_true = x$F, U_counts = U_counts, F_counts = F_counts, U_hat = U_hat, F_hat = F_hat)
  })

  posterior_draws <- reactive({
    obs <- observed_data()
    nsim <- suppressWarnings(as.integer(input$nsim)); if (is.na(nsim) || nsim <= 0) nsim <- 500L
    U_draw <- array(0, dim = c(obs$k, obs$k, nsim))
    F_draw <- array(0, dim = c(obs$k, obs$k, nsim))
    U_struct_zero <- obs$U_true == 0
    F_struct_zero <- obs$F_true == 0
    fix_zero <- isTRUE(input$fix_zero_structural)
    for (r in seq_len(nsim)) {
      for (j in seq_len(obs$k)) {
        if (fix_zero) {
          active <- !U_struct_zero[, j]
          if (any(active)) {
            y_col <- obs$U_counts[active, j]
            d_col <- max(obs$n[j] - sum(y_col), 0)
            s_col <- rdirichlet1(c(y_col, d_col) + 1)
            U_draw[active, j, r] <- s_col[seq_along(y_col)]
          }
        } else {
          y_col <- obs$U_counts[, j]
          d_col <- max(obs$n[j] - sum(y_col), 0)
          s_col <- rdirichlet1(c(y_col, d_col) + 1)
          U_draw[, j, r] <- s_col[seq_len(obs$k)]
        }
      }
      for (j in seq_len(obs$k)) {
        for (i in seq_len(obs$k)) {
          if (fix_zero && F_struct_zero[i, j]) {
            F_draw[i, j, r] <- 0
          } else {
            F_draw[i, j, r] <- rgamma(1, shape = obs$F_counts[i, j] + 1, rate = obs$n[j])
          }
        }
      }
    }
    list(U = U_draw, F = F_draw, k = obs$k, nsim = nsim)
  })

  density_data <- reactive({
    obs <- observed_data(); fix_zero <- isTRUE(input$fix_zero_structural)
    u_cells <- expand_grid(i = seq_len(obs$k), j = seq_len(obs$k)) %>%
      mutate(matrix = "U", point = obs$U_hat[cbind(i, j)], n = obs$n[j], structural_zero = obs$U_true[cbind(i, j)] == 0) %>%
      mutate(dens = pmap(list(point, n, structural_zero), function(point, n, structural_zero) {
        if (fix_zero && structural_zero) return(tibble(x = seq(0, 1, by = 0.01), d = NA_real_))
        calc_density_beta(point, n)
      })) %>% unnest(dens)
    f_cells <- expand_grid(i = seq_len(obs$k), j = seq_len(obs$k)) %>%
      mutate(matrix = "F", point = obs$F_hat[cbind(i, j)], n = obs$n[j], structural_zero = obs$F_true[cbind(i, j)] == 0) %>%
      mutate(dens = pmap(list(point, n, structural_zero), function(point, n, structural_zero) {
        if (fix_zero && structural_zero) return(tibble(x = seq(0, 12, length.out = 400), d = NA_real_))
        calc_density_gamma(point, n)
      })) %>% unnest(dens)
    bind_rows(u_cells, f_cells) %>% mutate(to = factor(paste0("to ", i), levels = paste0("to ", seq_len(obs$k))), from = factor(paste0("from ", j), levels = paste0("from ", seq_len(obs$k))))
  })

  cell_draws <- reactive({
    post <- posterior_draws()
    u_tbl <- as.data.frame.table(post$U, responseName = "value") %>% transmute(matrix = "U", i = as.integer(Var1), j = as.integer(Var2), draw = as.integer(Var3), value = value)
    f_tbl <- as.data.frame.table(post$F, responseName = "value") %>% transmute(matrix = "F", i = as.integer(Var1), j = as.integer(Var2), draw = as.integer(Var3), value = value)
    bind_rows(u_tbl, f_tbl)
  })

  cell_ci <- reactive({
    cell_draws() %>% group_by(matrix, i, j) %>% summarize(low95 = quantile(value, 0.025, na.rm = TRUE), upp95 = quantile(value, 0.975, na.rm = TRUE), .groups = "drop")
  })

  cell_post_med <- reactive({
    cell_draws() %>% group_by(matrix, i, j) %>% summarize(post_med = median(value, na.rm = TRUE), .groups = "drop")
  })

  derived_draws <- reactive({
    post <- posterior_draws()
    map_dfr(seq_len(post$nsim), function(r) {
      U <- post$U[, , r, drop = TRUE]
      F <- post$F[, , r, drop = TRUE]
      A <- U + F
      spec <- spectral_summary(A)
      repro_cols <- which(colSums(F) > 0)
      start_stage <- if (length(repro_cols) > 0) min(repro_cols) else 1L
      tibble(draw = r, quantity = quantity_levels, value = c(spec$lambda, life_expectancy_from_start(U, start_stage), generation_time_safe(U, F), spec$damping))
    }) %>% mutate(quantity = factor(quantity, levels = quantity_levels))
  })

  point_estimates <- reactive({
    obs <- observed_data(); A <- obs$U_hat + obs$F_hat
    spec <- spectral_summary(A)
    repro_cols <- which(colSums(obs$F_hat) > 0)
    start_stage <- if (length(repro_cols) > 0) min(repro_cols) else 1L
    tibble(quantity = factor(quantity_levels, levels = quantity_levels), point = c(spec$lambda, life_expectancy_from_start(obs$U_hat, start_stage), generation_time_safe(obs$U_hat, obs$F_hat), spec$damping))
  })

  derived_post_med <- reactive({
    derived_draws() %>% group_by(quantity) %>% summarize(post_med = median(value, na.rm = TRUE), .groups = "drop")
  })

  cell_labels <- function(which_mat = c("U", "F")) {
    which_mat <- match.arg(which_mat)
    obs <- observed_data()
    expand_grid(i = seq_len(obs$k), j = seq_len(obs$k)) %>%
      mutate(
        num = if (which_mat == "U") obs$U_counts[cbind(i, j)] else obs$F_counts[cbind(i, j)],
        den = obs$n[j],
        rate = if (which_mat == "U") obs$U_hat[cbind(i, j)] else obs$F_hat[cbind(i, j)],
        lab = paste0(num, "/", den, " (", fmt_val(rate), ")"),
        to = factor(paste0("to ", i), levels = paste0("to ", seq_len(obs$k))),
        from = factor(paste0("from ", j), levels = paste0("from ", seq_len(obs$k)))
      )
  }

  plot_matrix <- function(which_mat = c("U", "F")) {
    which_mat <- match.arg(which_mat)
    dd <- density_data() %>% filter(matrix == which_mat)
    pts <- dd %>% distinct(to, from, i, j, point)
    ci <- cell_ci() %>% filter(matrix == which_mat)
    dd2 <- dd %>% left_join(ci %>% select(i, j, low95, upp95), by = c("i", "j")) %>% mutate(d95 = if_else(x >= low95 & x <= upp95, d, NA_real_))
    dd_plot <- dd2 %>% filter(is.finite(d), d > 0)
    dd95_plot <- dd2 %>% filter(is.finite(d95), d95 > 0)
    nonempty <- dd2 %>% group_by(i, j) %>% summarize(has_density = any(is.finite(d) & d > 0), .groups = "drop") %>% filter(has_density)
    pts_plot <- pts %>% inner_join(nonempty, by = c("i", "j"))
    pts_post <- pts %>% left_join(cell_post_med() %>% filter(matrix == which_mat), by = c("i", "j")) %>% inner_join(nonempty, by = c("i", "j"))

    labs_df <- cell_labels(which_mat) %>%
      inner_join(
        dd_plot %>% group_by(i, j, to, from) %>% summarize(xmin = min(x), xmax = max(x), ymax = max(d), .groups = "drop") %>% mutate(x_text = xmin + 0.03 * (xmax - xmin), y_text = ymax * 0.93),
        by = c("i", "j", "to", "from")
      )

    col_fill <- if (which_mat == "U") u_col else f_col
    subtitle <- if (which_mat == "U") "Rows = to stage, columns = from stage" else "Recruitment in top row: from stage j to stage 1"

    ggplot(dd2, aes(x = x, y = d)) +
      geom_ribbon(data = dd_plot, fill = col_fill, alpha = 0.25, aes(ymin = 0, ymax = d)) +
      geom_ribbon(data = dd95_plot, fill = col_fill, alpha = 0.35, aes(ymin = 0, ymax = d95)) +
      geom_label(data = labs_df, aes(x = x_text, y = y_text, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.4, label.size = 0.15, fill = "white", alpha = 0.85) +
      geom_vline(data = pts_plot, aes(xintercept = point), linetype = 2, linewidth = 0.9, color = point_col) +
      geom_vline(data = pts_post, aes(xintercept = post_med), linetype = 3, linewidth = 0.6, color = post_col) +
      facet_grid(to ~ from, scales = "free_x", switch = "y") +
      labs(x = NULL, y = NULL, title = paste(which_mat, "matrix sampling distributions"), subtitle = subtitle) +
      theme_app_plot()
  }

  output$u_plot <- renderPlot({
    validate(need(validity()$ok, "Check inputs: values must be non-negative and sample sizes must be > 0."))
    plot_matrix("U")
  })

  output$f_plot <- renderPlot({
    validate(need(validity()$ok, "Check inputs: values must be non-negative and sample sizes must be > 0."))
    plot_matrix("F")
  })

  output$derived_plot <- renderPlot({
    validate(need(validity()$ok, "Check inputs: values must be non-negative and sample sizes must be > 0."))
    dd <- derived_draws(); pe <- point_estimates(); pm <- derived_post_med()
    ci <- dd %>% group_by(quantity) %>% summarize(low95 = quantile(value, 0.025, na.rm = TRUE), upp95 = quantile(value, 0.975, na.rm = TRUE), .groups = "drop")
    dens <- dd %>% group_by(quantity) %>% group_modify(~ {
      x <- .x$value[is.finite(.x$value)]
      if (length(x) < 10) return(tibble(x = NA_real_, d = NA_real_))
      kd <- density(x, n = 512, na.rm = TRUE)
      tibble(x = kd$x, d = kd$y)
    }) %>% ungroup() %>% left_join(ci, by = "quantity") %>% mutate(d95 = if_else(x >= low95 & x <= upp95, d, NA_real_))

    ggplot(dens, aes(x = x, y = d)) +
      geom_ribbon(data = dens %>% filter(is.finite(d), d > 0), aes(ymin = 0, ymax = d), fill = derived_col, alpha = 0.25) +
      geom_ribbon(data = dens %>% filter(is.finite(d95), d95 > 0), aes(ymin = 0, ymax = d95), fill = derived_col, alpha = 0.35) +
      geom_vline(data = pe, aes(xintercept = point), linetype = 2, linewidth = 0.95, color = point_col) +
      geom_vline(data = pm, aes(xintercept = post_med), linetype = 3, linewidth = 0.65, color = post_col) +
      facet_wrap(~ quantity, ncol = 1, scales = "free") +
      labs(x = NULL, y = "Density", title = "Derived quantities") +
      theme_app_plot()
  })

  output$summary_table <- renderTable({
    validate(need(validity()$ok, "Check inputs."))
    dd <- derived_draws(); pe <- point_estimates()
    dd %>% group_by(quantity) %>% summarize(med = median(value, na.rm = TRUE), low95 = quantile(value, 0.025, na.rm = TRUE), upp95 = quantile(value, 0.975, na.rm = TRUE), .groups = "drop") %>%
      left_join(pe, by = "quantity") %>% mutate(across(where(is.numeric), ~ signif(.x, 4))) %>%
      transmute(`Derived quantity` = quantity, `Point estimate` = point, `Posterior median` = med, `Lower 95% bound` = low95, `Upper 95% bound` = upp95)
  }, digits = 4)

  output$details_ui <- renderUI({
    if (!details_open()) return(NULL)
    obs <- tryCatch(observed_data(), error = function(e) NULL)
    setup_txt <- if (is.null(obs)) "Current setup unavailable." else paste0("Dimension k = ", obs$k, "; n = (", paste(obs$n, collapse = ", "), ")")

    tagList(
      div(class = "details-overlay"),
      div(class = "details-panel",
          div(style = "display:flex; justify-content:space-between; align-items:center;", tags$h4("Details"), actionButton("close_details", "Close")),
          tags$p("This app mirrors the manuscript workflow for uncertainty propagation in stage-structured matrix population models."),
          tags$h5("How to use"),
          tags$ul(
            tags$li("Pick a dimension tab (2x2, 3x3, or 4x4)."),
            tags$li("Enter stage sample sizes n[j] for each source stage."),
            tags$li("Enter observed transition counts in U and observed recruit counts in F."),
            tags$li("Use the matrix panels to inspect cell-level uncertainty and the right panel/table to inspect derived quantities.")
          ),
          tags$h5("Model setup"),
          tags$ul(
            tags$li("A = U + F, rows = destination stage, columns = source stage."),
            tags$li("U entries represent surviving transitions; deaths are implicit as the remaining count in each source stage."),
            tags$li("F entries represent recruit counts per source stage."),
            tags$li("U uncertainty uses Dirichlet sampling on living transitions plus death residual."),
            tags$li("F uncertainty uses a Gamma posterior for Poisson count-rate sampling."),
            tags$li("Structural-zero mode can fix entered zeros at zero in posterior draws."),
            tags$li("Dashed red lines mark point estimates; dotted grey lines mark posterior medians."),
            tags$li(paste("Posterior draws currently set to:", input$nsim))
          ),
          tags$p(tags$b("Current setup: "), setup_txt)
      )
    )
  })
}

shinyApp(ui, server)
