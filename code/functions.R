# Shared helper functions for analyses and plotting.

# Required libraries ######################################################### ----
# require(rstan)
# require(loo)
# require(tidyverse)



# Transformations ############################################################ ----
# Convert probabilities in (0,1) to log-odds.
logit <- function(p) {
  log(p / (1 - p))
}

# Convert log-odds to probabilities in (0,1).
logit_inv <- function(a) {
  exp(a) / (exp(a) + 1)
}

# Convert a numeric vector to simplex proportions.
softmax <- function(x) {
  exp(x) / sum(exp(x))
}



# Distributions ############################################################## ----
# Draw one Dirichlet random vector (alpha is the concentration vector).
rdirichlet <- function(alpha) {
  M <- length(alpha)
  x <- rgamma(M, alpha)
  return(x / sum(x))
}


# Helpers to convert MPMs from Ellis et al. 2012 ############################# ----
# Parse a matrix encoded as text into a numeric square matrix.
string_to_mat <- function(A) {
  A <- gsub(pattern = "\\[|\\]|\\;", "", A)
  A <- strsplit(x = A, split = " ")[[1]]
  mat <- matrix(as.numeric(A), nrow = sqrt(length(A)), byrow = TRUE)
  return(mat)
}

# Parse stage sample sizes encoded as text into an integer vector.
nx_to_vec <- function(x) {
  x <- gsub("\\[|\\]", "", x)
  x <- strsplit(x, " ")[[1]]
  x <- gsub("NA", NA_integer_, x)
  x <- as.integer(x)
  return(x)
}



# Rcompadre helpers ########################################################## ----
# Return a compact view of a CompadreDB object with selected metadata columns.
cdb_glimpse <- function(db, cols = NULL) {
  db <- tibble::as_tibble(db)
  db <- dplyr::rename(db, StartYear = MatrixStartYear, EndYear = MatrixEndYear)
  db[, c("SpeciesAuthor", "MatrixPopulation", "MatrixComposite",
         "MatrixTreatment", "StartYear", "EndYear", cols)]
}

# Row-bind multiple CompadreDB objects while preserving class and version.
cdb_bind_rows <- function(dbs) {
  if (methods::is(dbs[[1]], "CompadreDB")) {
    vers <- dbs[[1]]@version
  } else {
    vers <- load_compadre(corrected = TRUE)@version
  }
  dbs <- dplyr::bind_rows(lapply(dbs, tibble::as_tibble))
  new("CompadreDB",
      data = dbs,
      version = vers)
}

# Wrap a tibble with a `mat` list-column as a CompadreDB object.
tibble_to_cdb <- function(tbl, version) {
  new("CompadreDB", data = tibble::as_tibble(tbl), version = version)
}

# Build a CompadreMat by averaging matching slots across rows.
compadre_mat_mean <- function(mats) {
  mats <- as.list(mats)
  if (length(mats) == 0) {
    stop("mats must contain at least one CompadreMat", call. = FALSE)
  }
  if (length(mats) == 1) {
    return(mats[[1]])
  }

  out <- mats[[1]]
  out@matU <- mat_mean(lapply(mats, matU))
  out@matF <- mat_mean(lapply(mats, matF))
  out@matC <- mat_mean(lapply(mats, matC))
  out@matA <- out@matU + out@matF + out@matC
  out
}

# Build one population-level matrix per group from the best available composite.
population_matrices_from_available <- function(db,
                                               preferred = c("Mean", "Individual"),
                                               group_cols = "MatrixPopulation") {
  db_tbl <- tibble::as_tibble(db)

  if ("MatrixPopulation" %in% names(db_tbl)) {
    pop_vals <- unique(db_tbl$MatrixPopulation)
    has_specific <- any(!grepl(";", pop_vals, fixed = TRUE)) &&
      any(grepl(";", pop_vals, fixed = TRUE))
    if (has_specific) {
      db_tbl <- db_tbl %>% dplyr::filter(!grepl(";", MatrixPopulation, fixed = TRUE))
    }

    has_pooled_label <- "Pooled" %in% pop_vals && any(pop_vals != "Pooled")
    if (has_pooled_label) {
      db_tbl <- db_tbl %>% dplyr::filter(MatrixPopulation != "Pooled")
    }
  }

  split_groups <- split(
    db_tbl,
    interaction(db_tbl[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  )

  purrr::map_dfr(split_groups, function(df_group) {
    use <- NULL

    for (comp in preferred) {
      cand <- df_group[df_group$MatrixComposite == comp, , drop = FALSE]
      if (nrow(cand) > 0) {
        use <- cand
        break
      }
    }

    if (is.null(use)) {
      return(df_group[0, , drop = FALSE])
    }

    out <- use[1, , drop = FALSE]
    start_year <- suppressWarnings(as.integer(use$MatrixStartYear))
    end_year <- suppressWarnings(as.integer(use$MatrixEndYear))

    out$MatrixComposite <- preferred[[1]]
    out$MatrixStartYear <- if (all(is.na(start_year))) out$MatrixStartYear[[1]] else min(start_year, na.rm = TRUE)
    out$MatrixEndYear <- if (all(is.na(end_year))) out$MatrixEndYear[[1]] else max(end_year, na.rm = TRUE)
    out$mat <- list(compadre_mat_mean(use$mat))
    out
  })
}



# Other utilities ############################################################ ----
# Load the first object from an .RData file and return it.
rdata_load <- function(path) {
  env <- new.env()
  x <- load(path, env)[1]
  return(env[[x]])
}

# Load the first object from an .RData file and drop unused columns.
rdata_load2 <- function(path) {
  env <- new.env()
  x <- load(path, env)[1]
  out <- env[[x]]
  if ("MatrixID" %in% names(out)) {
    out$MatrixID <- as.character(out$MatrixID)
  }
  out$Altitude <- NULL
  out$MatrixStartYear <- NULL
  out$MatrixEndYear <- NULL
  return(out)
}

# Wrapper for life expectancy from Rage with default start stage behaviour.
life_expect <- function(matU, mixdist = NULL, start = NULL) {
  if (is.null(mixdist) && is.null(start)) {
    start <- 1L
  }
  Rage::life_expect_mean(matU, mixdist = mixdist, start = start)
}

# Proportion of cohort expected to be in reproductive stages at first maturity.
repro_prop_start <- function(matU, start, rep_stages) {
  Rage::mature_distrib(matU, start = start, repro_stages = rep_stages)
}

# Survivorship trajectory starting from mature-stage mixture.
lx_from_mature <- function(matU, rep_prop1, xmax = 1000, lx_crit = -1) {
  Rage::mpm_to_lx(matU, start = rep_prop1, xmax = xmax, lx_crit = lx_crit)
}

# Time to quasi-stationary distribution for a projection matrix.
qsd <- function(matU, rep_prop1, conv = 0.01, nmax = 1e5) {
  Rage::qsd_converge(matU, rep_prop1, conv = conv, N = nmax)
}

# Safe qsd wrapper that returns NA when convergence fails.
qsd_safe <- function(matU, rep_prop1, conv = 0.01, nmax = 1e5) {
  tryCatch(
    Rage::qsd_converge(matU, rep_prop1, conv = conv, N = nmax),
    error = function(e) NA_integer_
  )
}

# Net reproductive rate from U and F submatrices.
R0 <- function(matU, matF) {
  Rage::net_repro_rate(matU, matF)
}

# Script helpers ############################################################# ----
# Collapse multiple source strings into a single ';'-separated unique string.
collapse_fn <- function(x) {
  ifelse(all(is.na(x)),
         NA_character_,
         paste(unique(x[!is.na(x)]), collapse = "; "))
}

# Plot helpers ############################################################### ----
# Generate a viridis palette used across manuscript figures.
mpm_pal <- function(n = 5, option = "viridis", begin = 0.1, end = 0.85) {
  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    stop("Missing package: viridisLite. Install with install.packages(\"viridisLite\").",
         call. = FALSE)
  }
  viridisLite::viridis(n, option = option, begin = begin, end = end)
}

# Named color list used for consistent point/line/fill assignments.
mpm_colors <- function(option = "viridis") {
  pal <- mpm_pal(5, option = option)
  list(
    light = pal[1],
    fill = pal[2],
    mid = pal[3],
    accent = pal[4],
    dark = pal[5]
  )
}

# Common ggplot theme used across analyses.
theme_mpm <- function(base_size = 11.5) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "grey80", linewidth = 0.4, fill = NA),
      strip.background = ggplot2::element_rect(color = "grey80", fill = "grey90", linewidth = 0.4),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(face = "plain"),
      legend.key = ggplot2::element_blank()
    )
}

# Apply common plotting defaults (theme and geom defaults) globally.
set_mpm_plot_defaults <- function(base_size = 11.5) {
  ggplot2::theme_set(theme_mpm(base_size = base_size))
  ggplot2::update_geom_defaults("line", list(linewidth = 0.5))
  ggplot2::update_geom_defaults("point", list(size = 1.2))
  ggplot2::update_geom_defaults("linerange", list(linewidth = 0.4))
  invisible(NULL)
}

# Extract PRISM temperature and precipitation at species coordinates.
fetch_prism <- function(file_tmp, file_ppt, spp) {
  prism_tmp <- terra::rast(file_tmp)
  prism_ppt <- terra::rast(file_ppt)
  if (!is.data.frame(spp)) spp <- spp[[1]]
  spp_sf <- sf::st_as_sf(spp, coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  spp_sf <- sf::st_transform(spp_sf, terra::crs(prism_tmp))
  tmp_vals <- terra::extract(prism_tmp, terra::vect(spp_sf))[, 2]
  ppt_vals <- terra::extract(prism_ppt, terra::vect(spp_sf))[, 2]
  out <- sf::st_drop_geometry(spp_sf)
  out$tmp <- tmp_vals
  out$ppt <- ppt_vals
  return(tibble::as_tibble(out))
}



# Confirm that stage-specific sample sizes match transition rates ############ ----
# Check all U columns for compatibility with observed stage sample sizes.
check_freqs_mat <- function(matU, N, prec = 0.001) {
  dim <- nrow(matU)
  out <- character(dim)

  for (i in seq_len(dim)) {
    out[i] <- check_freqs(N[i], matU[, i], prec = prec)
  }

  return(data.frame(N, x = out))
}

# Check whether proportions are representable by integer counts at sample size n.
check_freqs <- function(n, x, prec = 0.001) {
  if (is.na(n) || n == 0) {
    return(NA)
  } else {
    y <- vector(mode = "numeric", length = length(x))

    for (i in seq_along(x)) {
      y[i] <- round(x[i] * n) / n
    }

    check <- all(abs(y - x) <= prec)
    if (check) {
      return("Pass")
    } else {
      y <- round(y, as.integer(log10(1 / prec)))
      return(paste(paste(round(y[y > 0], 6), collapse = "; ")))
    }
  }
}



# Sampling distributions for single mpm ###################################### ----
# Posterior density over transition probability with flat prior.
# Uses binomial for survival transitions and Poisson for fecundity entries.
dens_fn <- function(x, n, fec) {
  # x is number of successes
  # n is number of trials
  if (is.na(x)) {
    return(data.frame(p = 0, pp = NA_real_))
  } else {
    p <- seq(0, 1, 0.01)    # population probability

    if (fec) {
      l <- dpois(x, p * n)
      cn <- integrate(function(z) dpois(x, z * n), lower = 0, upper = 1)$value
    } else {
      l <- dbinom(x, n, p)    # likelihood
      cn <- integrate(function(z) dbinom(x, n, z), lower = 0, upper = 1)$value
    }
    pp <- l / cn # posterior probability (assuming flat prior)
    return(data.frame(p = p, pp = pp))
  }
}

# Simulate one U-column from observed stage size using Dirichlet posterior.
sim_stage_U <- function(x, vital_ind, n) {
  colsum <- sum(x)
  if (length(vital_ind) == 0) {
    out <- x
  } else {
    if (colsum > 1) x <- x / colsum
    mortality <- 1 - sum(x)
    rates <- c(x[vital_ind], mortality)
    frequencies <- rates * n
    rates_sim <- rdirichlet(frequencies + 1)
    out <- numeric(length(x))
    out[vital_ind] <- rates_sim[seq_along(vital_ind)]
  }
  return(out)
}

# Simulate one F-column from observed stage size using gamma-Poisson update.
sim_stage_F <- function(x, vital_ind, n) {
  if (length(vital_ind) == 0 || is.na(n)) {
    out <- x
  } else {
    rates <- x[vital_ind]
    y <- rates * n
    rates_sim <- sapply(y, function(y) rgamma(1, shape = 1 + y, rate = 0 + n))
    out <- numeric(length(x))
    out[vital_ind] <- rates_sim
  }
  return(out)
}

# Simulate full U matrix conditional on column sample sizes.
sim_U <- function(matU, posU, N) {
  if ("list" %in% class(matU)) matU <- matU[[1]]
  if ("list" %in% class(N)) N <- unlist(N)
  N <- as.numeric(N)
  if (length(N) == 0) N <- rep(NA_real_, ncol(matU))
  if (length(N) < ncol(matU)) N <- c(N, rep(NA_real_, ncol(matU) - length(N)))
  if (length(N) > ncol(matU)) N <- N[seq_len(ncol(matU))]

  simU <- matrix(0, nrow = nrow(matU), ncol = ncol(matU))

  for (i in seq_len(ncol(matU))) {
    if (is.na(N[i])) {
      simU[, i] <- matU[, i]
    } else if (N[i] == 0) {
      simU[, i] <- NA_real_
    } else {
      vital_ind <- which(posU[, i] > 0)
      simU[, i] <- sim_stage_U(matU[, i], vital_ind, N[i])
    }
  }
  return(simU)
}

# Simulate full F matrix conditional on column sample sizes.
sim_F <- function(matF, posF, N) {
  if ("list" %in% class(matF)) matF <- matF[[1]]
  if ("list" %in% class(N)) N <- unlist(N)
  N <- as.numeric(N)
  if (length(N) == 0) N <- rep(NA_real_, ncol(matF))
  if (length(N) < ncol(matF)) N <- c(N, rep(NA_real_, ncol(matF) - length(N)))
  if (length(N) > ncol(matF)) N <- N[seq_len(ncol(matF))]

  simF <- matrix(0, nrow = nrow(matF), ncol = ncol(matF))

  for (i in seq_len(ncol(matF))) {
    if (is.na(N[i])) {
      simF[, i] <- matF[, i]
    } else if (N[i] == 0) {
      simF[, i] <- NA_real_
    } else {
      vital_ind <- which(posF[, i] > 0)
      simF[, i] <- sim_stage_F(matF[, i], vital_ind, N[i])
    }
  }
  return(simF)
}

# Generate nsim sampled U matrices.
sim_U_wrapper <- function(matU, posU = matU > 0, N, nsim) {
  return(replicate(nsim, sim_U(matU, posU, N), simplify = FALSE))
}

# Generate nsim sampled F matrices.
sim_F_wrapper <- function(matF, posF = matF > 0, N, nsim) {
  return(replicate(nsim, sim_F(matF, posF, N), simplify = FALSE))
}



# Rstan helpers ############################################################## ----
# Standard control settings used for Stan model fitting.
ctrl1 <- list(adapt_delta = 0.95, stepsize = 0.05)
ctrl2 <- list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 11)
ctrl3 <- list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 12)


# Extract one parameter array from a stanfit object.
rstan_extract <- function(fit, var) {
  rstan::extract(fit, var)[[var]]
}

# Summarize key convergence diagnostics from a stanfit object.
stan_diagnostics <- function(fit) {
  n <- length(rstan_extract(fit, "lp__"))
  stan_summary <- as.data.frame(summary(fit)$summary)
  rhat_high <- length(which(stan_summary$Rhat > 1.1))
  n_eff_low <- length(which(stan_summary$n_eff / n < 0.1))
  mcse_high <- length(which(stan_summary$se_mean / stan_summary$sd > 0.1))
  n_diverg <- rstan::get_num_divergent(fit)
  return(tibble::tibble(rhat_high, n_eff_low, mcse_high, n_diverg))
}

# Fit Stan model with fallback to stricter controls if diagnostics are poor.
stanfn <- function(object, data, control = ctrl1, iter = 3000,
                   pars_excl = NULL, seed = 12345) {
  fit <- rstan::sampling(object = object, data = data, warmup = 2000,
                         iter = iter, thin = 2, chains = 2, control = control,
                         pars = pars_excl, include = FALSE, seed = seed)

  # if signs of poor convergence, re-fit with ctrl2
  if (any(stan_diagnostics(fit) > 0)) {
    fit <- rstan::sampling(object = object, data = data, warmup = 2000,
                           iter = iter, thin = 2, chains = 2, control = ctrl2,
                           pars = pars_excl, include = FALSE, seed = seed)
  }

  return(fit)
}

# Summarize posterior intervals for a vector-valued parameter.
posterior_vec <- function(fit, x, var, exp = FALSE) {
  var <- rstan::extract(fit, var)[[var]]
  fn <- ifelse(exp,
               function(x, q) exp(quantile(x, q)),
               function(x, q) quantile(x, q))
  return(
    tibble::tibble(
      x = x,
      med = apply(var, 2, fn, q = 0.500),
      low80 = apply(var, 2, fn, q = 0.10),
      upp80 = apply(var, 2, fn, q = 0.90),
      low95 = apply(var, 2, fn, q = 0.025),
      upp95 = apply(var, 2, fn, q = 0.975)
    )
  )
}

# Summarize posterior predictions yhat for plotting model fit.
summarize_yhat <- function(fit, label) {
  yhat <- rstan::extract(fit, "yhat")$yhat

  return(tibble::tibble(n = seq_len(ncol(yhat)),
                        y = y,
                        model = label,
                        yhat_med = apply(yhat, 2, function(x) quantile(x, 0.50)),
                        yhat_low90 = apply(yhat, 2, function(x) quantile(x, 0.05)),
                        yhat_upp90 = apply(yhat, 2, function(x) quantile(x, 0.95))))
}

# Summarize posterior regression coefficients across lags.
summarize_beta <- function(fit, label, wt = FALSE) {
  if (wt) {
    beta <- rstan::extract(fit, "beta_wt")$beta_wt
  } else {
    beta <- rstan::extract(fit, "beta")$beta
  }

  return(tibble::tibble(lag = seq_len(ncol(beta)),
                        model = label,
                        beta_med = apply(beta, 2, function(x) quantile(x, 0.50)),
                        beta_low95 = apply(beta, 2, function(x) quantile(x, 0.05)),
                        beta_upp95 = apply(beta, 2, function(x) quantile(x, 0.95))))
}

# Summarize model fit statistics (LOO/WAIC) and diagnostics.
summarize_fit <- function(fit, label) {
  ll <- extract_log_lik(fit)
  lppd <- sum(log(colMeans(exp(ll))))

  loo_mat <- suppressWarnings(loo(fit)$estimates)
  elpd_loo <- loo_mat[1, 1]
  elpd_loo_se <- loo_mat[1, 2]

  waic_mat <- suppressWarnings(waic(ll)$estimates)
  elpd_waic <- waic_mat[1, 1]
  elpd_waic_se <- waic_mat[1, 2]

  return(cbind(tibble::tibble(model = label),
               stan_diagnostics(fit),
               tibble::tibble(elpd_loo, elpd_loo_se, elpd_waic, elpd_waic_se)))
}

# Summarize out-of-sample predictions and test log-likelihood.
summarize_xval <- function(fit, label) {
  ll_test <- extract_log_lik(fit, "log_lik_test")
  lppd_test <- sum(log(colMeans(exp(ll_test))))

  yhat_test <- rstan::extract(fit, "yhat_test")$yhat_test
  yhat_test_median <- apply(yhat_test, 2, median)

  return(dplyr::bind_cols(tibble::tibble(model = label),
                          stan_diagnostics(fit),
                          tibble::tibble(lppd_test, yhat_test = yhat_test_median)))
}



# MPM manipulations ########################################################## ----
# Convert matrix entries to a long data frame with indices and stage labels.
mpm_flatten <- function(matA, matU, matF, matC, stage_names) {
  d <- nrow(matU)
  base_int <- expand.grid(to_col = seq_len(d), from_col = seq_len(d))
  base_name <- expand.grid(to_name = stage_names, from_name = stage_names)
  out <- cbind.data.frame(base_int, base_name)
  out$A <- c(matA)
  out$U <- c(matU)
  out$F <- c(matF)
  out$C <- c(matC)
  return(out)
}

# Rescale U columns with total survival > 1 back to sum 1.
scale_U <- function(matU) {
  out <- apply(matU, 2, function(x) {
    if (any(sum(x) > 1)) {
      x / sum(x)
    } else {
      x
    }
  })
  dimnames(out) <- dimnames(matU)
  return(out)
}

# Mean matrix across a list, optionally replacing NA with zero.
mat_mean2 <- function(l, na.rm = TRUE, replace_na = TRUE) {
  m <- list_mean(l, na.rm = na.rm)
  if (replace_na) m[is.na(m)] <- 0
  return(m)
}

# Mean matrix across a list of same-dimension matrices.
list_mean <- function(l, na.rm = TRUE) {
  arr <- simplify2array(l)
  m <- apply(arr, 1:2, function(x) mean(x, na.rm = na.rm))
  dimnames(m) <- dimnames(l[[1]])
  return(m)
}



# MPM and age-from-stage analyses ############################################ ----
# Truncate survivorship vector at tmax and optionally drop trailing zeros.
lx_submax <- function(lx, tmax, strip_zero = TRUE) {
  upp <- min(tmax, length(lx))
  lx <- lx[1L:upp] / lx[1L]
  if (strip_zero) lx <- lx[lx > 0]
  return(lx)
}

# Compute shape metric on truncated survivorship (safe for short vectors).
shape_surv2 <- function(lx, q) {
  upp <- min(q, length(lx))
  if (q < 4) {
    NA_real_
  } else {
    tryCatch(
      shape_surv(lx[1:upp]),
      error = function(e) NA_real_
    )
  }
}

# Sum with NA propagation if all values are missing.
sum2 <- function(x) {
  ifelse(all(is.na(x)), NA_real_, sum(x, na.rm = TRUE))
}

# Pool transition count vectors across replicate samples.
pool_counts <- function(nl) {
  if (is.null(nl) || length(nl) == 0) {
    return(NA_real_)
  }

  keep <- vapply(nl, function(x) !is.null(x) && length(x) > 0, logical(1))
  nl <- nl[keep]
  if (length(nl) == 0) {
    return(NA_real_)
  }

  X <- do.call(rbind, nl)
  if (is.null(dim(X))) {
    return(as.numeric(X))
  }

  return(apply(X, 2, sum2))
}

# Build transition matrix from row/col index table and value column name.
make_mat <- function(df, d, tr) {
  m <- matrix(0, nrow = d, ncol = d)
  m[cbind(df$row, df$col)] <- df[[tr]]
  return(m)
}






perturb_cust <- function(matU, matF, posU = matU > 0, posF = matF > 0,
                         exclude = NULL, type = "sensitivity") {

  # validate arguments
  type <- match.arg(type, c("sensitivity", "elasticity"))

  # matrix dimension
  m <- nrow(matU)

  # excluded stage classes
  posU[exclude, ] <- FALSE
  posU[, exclude] <- FALSE

  # combine components into matA
  matA <- matU + matF

  # lower and upper triangles (reflecting growth and retrogression)
  lwr <- upr <- matrix(FALSE, nrow = m, ncol = m)
  lwr[lower.tri(lwr)] <- TRUE
  upr[upper.tri(upr)] <- TRUE

  posStasi <- posU & diag(m)
  posRetro <- posU & upr
  posProgr <- posU & lwr

  if (type == "sensitivity") {

    pertMat <- popbio::sensitivity(matA)

    stasis <- ifelse(!any(posStasi), NA_real_, sum(pertMat[posStasi]))
    retro  <- ifelse(!any(posRetro), NA_real_, sum(pertMat[posRetro]))
    progr <- ifelse(!any(posProgr), NA_real_, sum(pertMat[posProgr]))
    fecund <- ifelse(!any(posF), NA_real_, sum(pertMat[posF]))

  } else {

    pertMat <- popbio::elasticity(matA)

    propU <- matU / matA
    propU[!posU] <- NA_real_
    propU[matA == 0 & posU] <- 1

    propProgr <- propRetro <- propU
    propProgr[upper.tri(propU, diag = TRUE)] <- NA
    propRetro[lower.tri(propU, diag = TRUE)] <- NA

    propStasi <- matrix(NA_real_, nrow = m, ncol = m)
    diag(propStasi) <- diag(propU)

    propF <- matF / matA
    propF[!posF] <- NA_real_
    propF[matA == 0 & posF] <- 1

    stasis <- sum_elast(pertMat, posStasi, propStasi)
    retro  <- sum_elast(pertMat, posRetro, propRetro)
    progr <- sum_elast(pertMat, posProgr, propProgr)
    fecund <- sum_elast(pertMat, posF, propF)
  }

  return(list(stasis = stasis,
              retro = retro,
              progr = progr,
              fecund = fecund))
}



# convenience function to sum elasticities given the perturbation matrix, the
#  matrix of possible transitions, and the matrix reflecting the proportional
#  contribution of the given process to the given element
sum_elast <- function(pert_mat, pos_mat, prop_mat) {
  ifelse(!any(pos_mat), NA_real_, sum(pert_mat * prop_mat, na.rm = TRUE))
}
