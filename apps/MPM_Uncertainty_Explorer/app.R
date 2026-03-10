# MPM Uncertainty Explorer (tabbed dimensions: 2x2, 3x3, 4x4, 5x5)

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
    F[1, 2] <- 10L
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
  if (k == 4) {
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
    F[1, 2] <- 6L
    F[1, 3] <- 11L
    F[1, 4] <- 16L
    return(list(n = n, U = U, F = F))
  }
  # k == 5
  n <- c(20L, 18L, 15L, 10L, 8L)
  U <- matrix(0L, 5, 5)
  U[1, 1] <- 6L
  U[2, 1] <- 8L
  U[2, 2] <- 7L
  U[3, 2] <- 6L
  U[3, 3] <- 6L
  U[4, 3] <- 4L
  U[4, 4] <- 5L
  U[5, 4] <- 3L
  U[5, 5] <- 5L
  F <- matrix(0L, 5, 5)
  F[1, 2] <- 5L
  F[1, 3] <- 9L
  F[1, 4] <- 13L
  F[1, 5] <- 18L
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

calc_density_gamma <- function(rate_hat, n, x_max = NULL) {
  if (length(n) != 1 || is.na(n) || !is.finite(n) || n <= 0) {
    return(tibble(x = seq(0, 12, length.out = 400), d = NA_real_))
  }
  if (length(rate_hat) != 1 || is.na(rate_hat) || !is.finite(rate_hat) || rate_hat < 0) {
    return(tibble(x = seq(0, 12, length.out = 400), d = NA_real_))
  }
  y_obs <- round(rate_hat * n)
  if (is.null(x_max) || !is.finite(x_max) || x_max <= 0) {
    x_max <- max(12, qgamma(0.995, shape = y_obs + 1, rate = n))
  }
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
id_a <- function(prefix, i, j) paste0(prefix, "_a_", i, "_", j)

is_f_cell <- function(i, j) {
  i == 1 && j > 1
}

panel_inputs <- function(prefix, k, defs) {
  A <- defs$U + defs$F
  tagList(
    tags$h4("Stage sample sizes"),
    tags$small("Number of individuals observed in each stage (column j)."),
    fluidRow(lapply(seq_len(k), function(j) {
      column(width = max(2, floor(12 / k)), numericInput(id_n(prefix, j), paste0("n[", j, "]"), value = defs$n[j], min = 1, step = 1))
    })),
    tags$h4("A matrix input"),
    tags$small("Enter observed counts in A. Top row (except [1,1]) is treated as F (pink). All other cells are treated as U (green)."),
    fluidRow(lapply(seq_len(k), function(j) column(width = max(2, floor(12 / k)), tags$b(paste("from", j))))),
    lapply(seq_len(k), function(i) {
      fluidRow(
        column(12, tags$small(tags$b(paste("to", i)))),
        lapply(seq_len(k), function(j) {
          cell_class <- if (is_f_cell(i, j)) "cell-f" else "cell-u"
          cell_tag <- if (is_f_cell(i, j)) tags$span(style = "color:#b03a86;", "F") else tags$span(style = "color:#1e7f56;", "U")
          column(
            width = max(2, floor(12 / k)),
            div(
              class = cell_class,
              numericInput(
                id_a(prefix, i, j),
                tagList("a[", i, ",", j, "] ", cell_tag),
                value = A[i, j],
                min = 0,
                step = 1
              )
            )
          )
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
      .cell-u .form-control { border-left: 4px solid #35b779; }
      .cell-f .form-control { border-left: 4px solid #e78ac3; }
      .details-overlay { position: fixed; inset: 0; background: rgba(15, 25, 20, 0.22); z-index: 2500; }
      .details-panel { position: fixed; top: 0; right: 0; width: 430px; max-width: 94vw; height: 100vh; overflow-y: auto; background: #ffffff; border-left: 1px solid #d7e3dd; box-shadow: -10px 0 30px rgba(31, 42, 36, 0.18); z-index: 2600; padding: 14px 16px 18px 16px; }
      .details-fab-wrap { position: fixed; right: 18px; bottom: 18px; z-index: 2400; }
      .details-fab-wrap .btn { border-radius: 999px; padding: 10px 16px; font-weight: 600; box-shadow: 0 8px 20px rgba(31, 42, 36, 0.2); }
    "))),
  titlePanel("MPM Uncertainty Explorer"),
  sidebarLayout(
    sidebarPanel(
      tags$p("Enter observed counts from a stage-structured population study, then explore how sampling uncertainty in the matrix entries propagates to key demographic quantities (λ, life expectancy, generation time, damping ratio)."),
      tags$p(tags$b("Step 1:"), " Select a matrix size tab below."),
      tags$p(tags$b("Step 2:"), " Enter stage sample sizes and observed counts in the A matrix."),
      tags$p(tags$b("Step 3:"), " Read the plots and summary table to compare point estimates with posterior uncertainty."),
      tabsetPanel(
        id = "dim_tab",
        selected = "3x3",
        tabPanel("2x2", panel_inputs("d2", 2, build_defaults(2))),
        tabPanel("3x3", panel_inputs("d3", 3, build_defaults(3))),
        tabPanel("4x4", panel_inputs("d4", 4, build_defaults(4))),
        tabPanel("5x5", panel_inputs("d5", 5, build_defaults(5)))
      ),
      checkboxInput("fix_zero_structural", "Fix entered zeros as structural zeros (held at 0 in all posterior draws)", value = TRUE),
      selectInput("nsim", "Number of posterior draws", choices = c("300", "500", "1000"), selected = "500"),
      actionButton("reset_defaults", "Reset active tab to defaults")
    ),
    mainPanel(
      fluidRow(
        column(8,
               div(class = "plot-card", plotOutput("a_plot", height = "760px"))
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
           "5x5" = list(prefix = "d5", k = 5L, defs = build_defaults(5)),
           list(prefix = "d3", k = 3L, defs = build_defaults(3)))
  })
  
  observeEvent(input$reset_defaults, {
    cfg <- active_cfg()
    for (j in seq_len(cfg$k)) updateNumericInput(session, id_n(cfg$prefix, j), value = cfg$defs$n[j])
    A_def <- cfg$defs$U + cfg$defs$F
    for (i in seq_len(cfg$k)) {
      for (j in seq_len(cfg$k)) {
        updateNumericInput(session, id_a(cfg$prefix, i, j), value = A_def[i, j])
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
    
    A <- matrix(0, nrow = cfg$k, ncol = cfg$k)
    for (i in seq_len(cfg$k)) {
      for (j in seq_len(cfg$k)) {
        aij <- suppressWarnings(as.numeric(input[[id_a(cfg$prefix, i, j)]]))
        def_aij <- cfg$defs$U[i, j] + cfg$defs$F[i, j]
        A[i, j] <- if (is.null(aij) || !is.finite(aij)) def_aij else aij
      }
    }
    list(k = cfg$k, n = n_vec, A = A)
  })
  
  validity <- reactive({
    x <- inputs()
    a_ok <- all(is.finite(x$A)) && all(x$A >= 0)
    n_ok <- all(is.finite(x$n)) && all(x$n > 0)
    list(ok = a_ok && n_ok)
  })
  
  observed_data <- reactive({
    req(validity()$ok)
    x <- inputs()
    is_f <- outer(seq_len(x$k), seq_len(x$k), Vectorize(is_f_cell))
    is_u <- !is_f
    A_counts <- matrix(pmax(as.integer(round(x$A)), 0L), nrow = x$k, ncol = x$k)

    U_raw <- A_counts
    U_raw[is_f] <- 0L
    F_raw <- matrix(0L, nrow = x$k, ncol = x$k)
    F_raw[is_f] <- A_counts[is_f]

    U_counts <- matrix(0, nrow = x$k, ncol = x$k)
    F_counts <- matrix(0, nrow = x$k, ncol = x$k)
    for (j in seq_len(x$k)) U_counts[, j] <- snap_u_count_column(U_raw[, j], x$n[j])
    F_counts[is_f] <- F_raw[is_f]

    U_hat <- sweep(U_counts, 2, x$n, "/")
    F_hat <- sweep(F_counts, 2, x$n, "/")
    A_hat <- U_hat + F_hat

    list(
      k = x$k, n = x$n, is_f = is_f, is_u = is_u,
      A_true = x$A, U_true = U_raw, F_true = F_raw,
      U_counts = U_counts, F_counts = F_counts,
      U_hat = U_hat, F_hat = F_hat, A_hat = A_hat,
      A_counts = A_counts
    )
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
          active <- obs$is_u[, j] & !U_struct_zero[, j]
          if (any(active)) {
            y_col <- obs$U_counts[active, j]
            d_col <- max(obs$n[j] - sum(y_col), 0)
            s_col <- rdirichlet1(c(y_col, d_col) + 1)
            U_draw[active, j, r] <- s_col[seq_along(y_col)]
          }
        } else {
          active <- obs$is_u[, j]
          y_col <- obs$U_counts[active, j]
          d_col <- max(obs$n[j] - sum(y_col), 0)
          s_col <- rdirichlet1(c(y_col, d_col) + 1)
          U_draw[active, j, r] <- s_col[seq_along(y_col)]
        }
      }
      for (j in seq_len(obs$k)) {
        for (i in seq_len(obs$k)) {
          if (!obs$is_f[i, j]) {
            F_draw[i, j, r] <- 0
          } else if (fix_zero && F_struct_zero[i, j]) {
            F_draw[i, j, r] <- 0
          } else {
            F_draw[i, j, r] <- rgamma(1, shape = obs$F_counts[i, j] + 1, rate = obs$n[j])
          }
        }
      }
    }
    list(U = U_draw, F = F_draw, A = U_draw + F_draw, k = obs$k, nsim = nsim)
  })
  
  density_data <- reactive({
    obs <- observed_data(); fix_zero <- isTRUE(input$fix_zero_structural)
    f_cells_meta <- expand_grid(i = seq_len(obs$k), j = seq_len(obs$k)) %>%
      filter(obs$is_f[cbind(i, j)]) %>%
      mutate(
        point = obs$A_hat[cbind(i, j)],
        n = obs$n[j],
        y_obs = round(point * n),
        q995 = qgamma(0.995, shape = y_obs + 1, rate = n)
      )
    f_xmax_common <- if (nrow(f_cells_meta) > 0) max(12, max(f_cells_meta$q995, na.rm = TRUE)) else 12

    expand_grid(i = seq_len(obs$k), j = seq_len(obs$k)) %>%
      mutate(
        cell_type = if_else(obs$is_f[cbind(i, j)], "F", "U"),
        point = obs$A_hat[cbind(i, j)],
        n = obs$n[j],
        structural_zero = if_else(cell_type == "F", obs$F_true[cbind(i, j)] == 0, obs$U_true[cbind(i, j)] == 0)
      ) %>%
      mutate(dens = pmap(list(point, n, structural_zero, cell_type), function(point, n, structural_zero, cell_type) {
        if (fix_zero && structural_zero) {
          if (cell_type == "F") return(tibble(x = seq(0, f_xmax_common, length.out = 400), d = NA_real_))
          return(tibble(x = seq(0, 1, by = 0.01), d = NA_real_))
        }
        if (cell_type == "F") return(calc_density_gamma(point, n, x_max = f_xmax_common))
        calc_density_beta(point, n)
      })) %>%
      unnest(dens) %>%
      mutate(
        to = factor(paste0("to ", i), levels = paste0("to ", seq_len(obs$k))),
        from = factor(paste0("from ", j), levels = paste0("from ", seq_len(obs$k)))
      )
  })
  
  cell_draws <- reactive({
    post <- posterior_draws()
    as.data.frame.table(post$A, responseName = "value") %>%
      transmute(i = as.integer(Var1), j = as.integer(Var2), draw = as.integer(Var3), value = value)
  })
  
  cell_ci <- reactive({
    cell_draws() %>% group_by(i, j) %>% summarize(low95 = quantile(value, 0.025, na.rm = TRUE), upp95 = quantile(value, 0.975, na.rm = TRUE), .groups = "drop")
  })
  
  cell_post_med <- reactive({
    cell_draws() %>% group_by(i, j) %>% summarize(post_med = median(value, na.rm = TRUE), .groups = "drop")
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
  
  cell_labels <- function() {
    obs <- observed_data()
    expand_grid(i = seq_len(obs$k), j = seq_len(obs$k)) %>%
      mutate(
        num = obs$A_counts[cbind(i, j)],
        den = obs$n[j],
        rate = obs$A_hat[cbind(i, j)],
        lab = paste0(num, "/", den, " (", fmt_val(rate), ")"),
        cell_type = if_else(obs$is_f[cbind(i, j)], "F", "U"),
        to = factor(paste0("to ", i), levels = paste0("to ", seq_len(obs$k))),
        from = factor(paste0("from ", j), levels = paste0("from ", seq_len(obs$k)))
      )
  }
  
  plot_matrix_a <- function() {
    dd <- density_data()
    obs <- observed_data()
    panel_levels <- unlist(lapply(seq_len(obs$k), function(i) {
      paste0("to ", i, " | from ", seq_len(obs$k))
    }))
    pts <- dd %>% distinct(to, from, i, j, point, cell_type)
    ci <- cell_ci()
    dd2 <- dd %>% left_join(ci %>% select(i, j, low95, upp95), by = c("i", "j")) %>% mutate(d95 = if_else(x >= low95 & x <= upp95, d, NA_real_))
    dd2 <- dd2 %>% mutate(panel = factor(paste0("to ", i, " | from ", j), levels = panel_levels))
    dd_plot <- dd2 %>% filter(is.finite(d), d > 0)
    dd95_plot <- dd2 %>% filter(is.finite(d95), d95 > 0)
    nonempty <- dd2 %>% group_by(i, j) %>% summarize(has_density = any(is.finite(d) & d > 0), .groups = "drop") %>% filter(has_density)
    pts_plot <- pts %>% inner_join(nonempty, by = c("i", "j")) %>% mutate(panel = factor(paste0("to ", i, " | from ", j), levels = panel_levels))
    pts_post <- pts %>% left_join(cell_post_med(), by = c("i", "j")) %>% inner_join(nonempty, by = c("i", "j")) %>% mutate(panel = factor(paste0("to ", i, " | from ", j), levels = panel_levels))
    
    labs_df <- cell_labels() %>%
      inner_join(
        dd_plot %>% group_by(i, j, to, from) %>% summarize(xmin = min(x), xmax = max(x), ymax = max(d), .groups = "drop") %>% mutate(x_text = xmin + 0.03 * (xmax - xmin), y_text = ymax * 0.93),
        by = c("i", "j", "to", "from")
      ) %>%
      mutate(panel = factor(paste0("to ", i, " | from ", j), levels = panel_levels))
    
    panel_limits <- expand_grid(i = seq_len(obs$k), j = seq_len(obs$k)) %>%
      mutate(
        cell_type = if_else(obs$is_f[cbind(i, j)], "F", "U"),
        xmin = 0,
        xmax = if_else(cell_type == "U", 1, max(dd2$x[dd2$cell_type == "F"], na.rm = TRUE)),
        to = factor(paste0("to ", i), levels = paste0("to ", seq_len(obs$k))),
        from = factor(paste0("from ", j), levels = paste0("from ", seq_len(obs$k))),
        panel = factor(paste0("to ", i, " | from ", j), levels = panel_levels)
      )

    ggplot(dd2, aes(x = x, y = d)) +
      geom_blank(data = panel_limits, aes(x = xmin, y = 0), inherit.aes = FALSE) +
      geom_blank(data = panel_limits, aes(x = xmax, y = 0), inherit.aes = FALSE) +
      geom_ribbon(data = dd_plot, aes(ymin = 0, ymax = d, fill = cell_type), alpha = 0.25) +
      geom_ribbon(data = dd95_plot, aes(ymin = 0, ymax = d95, fill = cell_type), alpha = 0.35) +
      geom_label(data = labs_df, aes(x = x_text, y = y_text, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.4, label.size = 0.15, fill = "white", alpha = 0.85) +
      geom_vline(data = pts_plot, aes(xintercept = point), linetype = 2, linewidth = 0.9, color = point_col) +
      geom_vline(data = pts_post, aes(xintercept = post_med), linetype = 3, linewidth = 0.6, color = post_col) +
      facet_grid(to ~ from, scales = "free_x", switch = "y") +
      scale_fill_manual(values = c("U" = u_col, "F" = f_col), breaks = c("U", "F"), labels = c("U matrix entry", "F matrix entry")) +
      labs(
        x = NULL, y = NULL,
        title = "A matrix — sampling distributions",
        subtitle = "Top row except [1,1] treated as F; all other cells treated as U",
        fill = "Cell type"
      ) +
      facet_wrap(~ panel, ncol = obs$k, scales = "free_x") +
      theme_app_plot() +
      theme(
        legend.position = "bottom",
        strip.text = element_text(size = 8),
        panel.spacing = grid::unit(0.7, "lines")
      )
  }
  
  output$a_plot <- renderPlot({
    validate(need(validity()$ok, "Check inputs: all values must be non-negative and sample sizes must be > 0."))
    plot_matrix_a()
  })
  
  output$derived_plot <- renderPlot({
    validate(need(validity()$ok, "Check inputs: all values must be non-negative and sample sizes must be > 0."))
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
      left_join(pe, by = "quantity") %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
      transmute(`Derived quantity` = quantity, `Point estimate` = point, `Posterior median` = med, `Lower 95% bound` = low95, `Upper 95% bound` = upp95)
  }, digits = 2)
  
  output$details_ui <- renderUI({
    if (!details_open()) return(NULL)
    obs <- tryCatch(observed_data(), error = function(e) NULL)
    setup_txt <- if (is.null(obs)) "Current setup unavailable." else paste0("k = ", obs$k, "; n = (", paste(obs$n, collapse = ", "), ")")
    
    tagList(
      div(class = "details-overlay"),
      div(class = "details-panel",
          div(style = "display:flex; justify-content:space-between; align-items:center;", tags$h4("About this app"), actionButton("close_details", "Close")),
          
          tags$p("This app illustrates how sampling uncertainty in matrix population model (MPM) entries propagates to key demographic outputs. Most published analyses treat matrix entries as exact — here, each entry is estimated from finite counts, and that uncertainty flows through to quantities like λ and generation time."),
          
          tags$h5("Workflow"),
          tags$ol(
            tags$li("Choose a matrix dimension tab (2×2, 3×3, 4×4, or 5×5)."),
            tags$li(HTML("Enter <b>n[j]</b>: the number of individuals observed in each source stage j.")),
            tags$li(HTML("Enter observed counts in the <b>A matrix</b>.")),
            tags$li("The app maps A to submatrices internally: top row except [1,1] is treated as F, and all other cells are treated as U."),
            tags$li("Adjust the posterior draws setting and the structural-zero option as needed, then read the plots.")
          ),
          
          tags$h5("Statistical model"),
          tags$ul(
            tags$li(HTML("<b>Projection matrix:</b> A = U + F, where rows index destination stage and columns index source stage.")),
            tags$li(HTML("<b>U columns:</b> Sampled from a Dirichlet posterior over (survive to stage 1, survive to stage 2, …, die). Deaths are implicit — the residual count in each column after accounting for all observed transitions.")),
            tags$li(HTML("<b>F cells:</b> Sampled from a Gamma posterior, treating recruit counts as Poisson with unknown rate. The posterior is Gamma(y + 1, n[j]), where y is the observed count.")),
            tags$li(HTML("<b>Structural zeros:</b> When enabled, cells entered as zero are fixed at zero in every draw — they represent transitions that are biologically impossible, not just unobserved."))
          ),
          
          tags$h5("Reading the plots"),
          tags$ul(
            tags$li(HTML("<span style='color:#c7352a;'>&#8211;&#8211; Red dashed line:</span> point estimate (observed rate).")),
            tags$li(HTML("<span style='color:#7a7a7a;'>&#183;&#183;&#183; Grey dotted line:</span> posterior median.")),
            tags$li("Darker shading shows the central 95% posterior interval; lighter shading shows the full distribution."),
            tags$li("A wide distribution means high uncertainty — typically from a small sample size. Narrow distributions mean the estimate is well-constrained."),
            tags$li("If the posterior median diverges from the point estimate, finite-sample bias is present: the posterior accounts for asymmetry and boundary effects that the raw rate ignores.")
          ),
          
          tags$h5("Derived quantities"),
          tags$ul(
            tags$li(HTML("<b>Lambda (λ):</b> Asymptotic population growth rate — the dominant eigenvalue of A.")),
            tags$li(HTML("<b>Mature life expectancy (L):</b> Expected total time spent alive, starting from the first reproductive stage.")),
            tags$li(HTML("<b>Generation time (T):</b> Mean age of parents of newborns in a stable population (requires the Rage package).")),
            tags$li(HTML("<b>Damping ratio (ρ):</b> Ratio of the two largest eigenvalue moduli — higher values mean faster convergence to stable stage structure."))
          ),
          
          tags$p(tags$b("Active setup: "), setup_txt),
          tags$p(tags$b("Posterior draws: "), input$nsim)
      )
    )
  })
}

shinyApp(ui, server)
