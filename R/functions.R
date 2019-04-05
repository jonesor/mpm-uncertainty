

### Required libraries #########################################################
# require(rstan)
# require(loo)
# require(tidyverse)



### Transformations ############################################################
logit <- function(p) {
  log(p / (1 - p))
}

logit_inv <- function(a) {
  exp(a) / (exp(a) + 1)
}

softmax <- function(x) {
  exp(x) / sum(exp(x))
}



### Distributions ##############################################################
rdirichlet <- function(alpha) {
  M <- length(alpha)
  x <- rgamma(M, alpha)
  return(x / sum(x))
}



### Helpers to convert MPMs from Ellis et al. 2012 #############################
stringToMat <- function(A) {
  A <- gsub(pattern = "\\[|\\]|\\;", "", A)
  A <- strsplit(x = A, split = " ")[[1]]
  mat <- matrix(as.numeric(A), nrow = sqrt(length(A)), byrow = TRUE)
  return(mat)
}


NxToVec <- function(x) {
  x <- gsub('\\[|\\]', '', x)
  x <- strsplit(x, ' ')[[1]]
  x <- gsub('NA', NA_integer_, x)
  x <- as.integer(x)
  return(x)
}



### Rcompadre helpers ##########################################################
cdb_glimpse <- function(db, cols = NULL) {
  db <- as_tibble(db)
  db <- rename(db, StartYear = MatrixStartYear, EndYear = MatrixEndYear)
  db[,c("SpeciesAuthor", "MatrixPopulation", "MatrixComposite",
        "MatrixTreatment", "StartYear", "EndYear", cols)]
}

cdb_bind_rows <- function(dbs) {
  vers <- dbs[[1]]@version
  dbs <- bind_rows(lapply(dbs, as_tibble))
  new("CompadreDB",
      data = dbs,
      version = vers)
}



### Other utilities ############################################################
rdata_load <- function(path) {
  env <- new.env()
  x <- load(path, env)[1]
  return(env[[x]])
}



### Confirm that stage-specific sample sizes match transition rates ############
CheckFreqsMat <- function(matU, N, prec = 0.001) {
  dim <- nrow(matU)
  out <- vector(mode = 'character', length = dim)
  
  for (i in 1:dim) {
    out[i] <- CheckFreqs(N[i], matU[,i], prec = prec) 
  }
  
  return(data.frame(N, x = out))
}

CheckFreqs <- function(n, x, prec = 0.001) {
  if(is.na(n) | n == 0) {
    return(NA)
  } else {
    y <- vector(mode = 'numeric', length = length(x))
    
    for(i in 1:length(x)) {
      y[i] <- round(x[i] * n) / n
    }
    
    check <- all(abs(y - x) <= prec)
    if(check == TRUE) {
      return(paste('Pass'))
    } else {
      y <- round(y, as.integer(log10(1/prec)))
      return(paste(paste(round(y[y > 0], 6), collapse = '; ')))
    }
  }
}



### Sampling distributions for single mpm ######################################
dens_fn <- function(x, n, fec) {
  # x is number of successes
  # n is number of trials
  if(is.na(x)) return(data.frame(p = 0, pp = NA_real_))
  p <- seq(0, 1, 0.01)    # population probability
  
  if (fec == TRUE) {
    l <- dpois(x, p * n)
    cn <- integrate(function (z) dpois(x, z * n), lower = 0, upper = 1)$value
  } else {
    l <- dbinom(x, n, p)    # likelihood
    cn <- integrate(function (z) dbinom(x, n, z), lower = 0, upper = 1)$value
  }
  pp <- l / cn # posterior probability (assuming flat prior)
  return(data.frame(p = p, pp = pp))
}


SimStageMatU <- function(x, vital_ind, n) {
  colsum <- sum(x)
  if(length(vital_ind) == 0) {
    out <- x
  } else {
    if(colsum > 1) x <- x / colsum
    mortality <- 1 - sum(x)
    rates <- c(x[vital_ind], mortality)
    frequencies <- rates * n
    rates_sim <- rdirichlet(frequencies + 1)
    out <- numeric(length(x))
    out[vital_ind] <- rates_sim[1:length(vital_ind)]
  }
  return(out)
}


SimStageMatF <- function(x, vital_ind, n) {
  if(length(vital_ind) == 0 | is.na(n)) {
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


SimMatU <- function(matU, posU, N) {
  if (class(matU) == "list") matU <- matU[[1]]
  if (class(N) == "list") N <- unlist(N)
  
  simU <- matrix(0, nrow = nrow(matU), ncol = ncol(matU))
  
  for (i in 1:ncol(matU)) {
    if (is.na(N[i])) {
      simU[,i] <- matU[,i]
    } else if (N[i] == 0) {
      simU[,i] <- NA_real_
    } else {
      vital_ind <- which(posU[,i] > 0)
      simU[,i] <- SimStageMatU(matU[,i], vital_ind, N[i])
    }
  }
  return(simU)
}


SimMatF <- function(matF, posF, N) {
  if (class(matF) == "list") matF <- matF[[1]]
  if (class(N) == "list") N <- unlist(N)
  
  simF <- matrix(0, nrow = nrow(matF), ncol = ncol(matF))
  
  for (i in 1:ncol(matF)) {
    if (is.na(N[i])) {
      simF[,i] <- matF[,i]
    } else if (N[i] == 0) {
      simF[,i] <- NA_real_
    } else {
      vital_ind <- which(posF[,i] > 0)
      simF[,i] <- SimStageMatF(matF[,i], vital_ind, N[i])
    }
  }
  return(simF)
}


SimMatUWrapper <- function(matU, posU = matU > 0, N, nsim) {
  return(replicate(nsim, SimMatU(matU, posU, N), simplify = FALSE))
}


SimMatFWrapper <- function(matF, posF = matF > 0, N, nsim) {
  return(replicate(nsim, SimMatF(matF, posF, N), simplify = FALSE))
}



### Rstan helpers ##############################################################
ctrl1 <- list(adapt_delta = 0.95, stepsize = 0.05)
ctrl2 <- list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 11)
ctrl3 <- list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 12)


rstan_extract <- function(fit, var) {
  rstan::extract(fit, var)[[var]]
}


stan_diagnostics <- function(fit) {
  n <- length(rstan_extract(fit, "lp__"))
  stan_summary <- as.data.frame(summary(fit)$summary)
  rhat_high <- length(which(stan_summary$Rhat > 1.1))
  n_eff_low <- length(which(stan_summary$n_eff / n < 0.1))
  mcse_high <- length(which(stan_summary$se_mean / stan_summary$sd > 0.1))
  n_diverg <- rstan::get_num_divergent(fit)
  return(tibble(rhat_high, n_eff_low, mcse_high, n_diverg))
}


stanfn <- function(object, data, control = ctrl1, iter = 3000,
                   pars_excl = NULL) {
  fit <- sampling(object = object, data = data, warmup = 2000,
                  iter = iter, thin = 2, chains = 2, control = control,
                  pars = pars_excl, include = FALSE)
  
  # if signs of poor convergence, re-fit with ctrl2
  if (any(stan_diagnostics(fit) > 0)) {
    fit <- sampling(object = object, data = data, warmup = 2000,
                    iter = iter, thin = 2, chains = 2, control = ctrl2,
                    pars = pars_excl, include = FALSE)
  }
  
  return(fit)
}


posterior_vec <- function(fit, x, var, exp = FALSE) {
  var <- rstan::extract(fit, var)[[var]]
  fn <- ifelse(exp,
               function(x, q) exp(quantile(x, q)),
               function(x, q) quantile(x, q))
  tibble(
    x = x,
    med = apply(var, 2, fn, q = 0.500),
    low80 = apply(var, 2, fn, q = 0.10),
    upp80 = apply(var, 2, fn, q = 0.90),
    low95 = apply(var, 2, fn, q = 0.025),
    upp95 = apply(var, 2, fn, q = 0.975))
}


summarize_yhat <- function(fit, label) {
  yhat <- rstan::extract(fit, 'yhat')$yhat
  
  tibble(n = seq_len(ncol(yhat)),
         y = y,
         model = label,
         yhat_med = apply(yhat, 2, function(x) quantile(x, 0.50)),
         yhat_low90 = apply(yhat, 2, function(x) quantile(x, 0.05)),
         yhat_upp90 = apply(yhat, 2, function(x) quantile(x, 0.95)))
}


summarize_beta <- function(fit, label, wt = FALSE) {
  if (wt == TRUE) {
    beta <- rstan::extract(fit, 'beta_wt')$beta_wt
  } else {
    beta <- rstan::extract(fit, 'beta')$beta
  }
  
  tibble(lag = seq_len(ncol(beta)),
         model = label,
         beta_med = apply(beta, 2, function(x) quantile(x, 0.50)),
         beta_low95 = apply(beta, 2, function(x) quantile(x, 0.05)),
         beta_upp95 = apply(beta, 2, function(x) quantile(x, 0.95)))
}


summarize_fit <- function(fit, label) {
  ll <- extract_log_lik(fit)
  lppd <- sum(log(colMeans(exp(ll))))
  
  loo_mat <- suppressWarnings(loo(fit)$estimates)
  elpd_loo <- loo_mat[1,1]
  elpd_loo_se <- loo_mat[1,2]
  
  waic_mat <- suppressWarnings(waic(ll)$estimates)
  elpd_waic <- waic_mat[1,1]
  elpd_waic_se <- waic_mat[1,2]
  
  cbind(tibble(model = label),
        stan_diagnostics(fit),
        tibble(elpd_loo, elpd_loo_se, elpd_waic, elpd_waic_se))
}


summarize_xval <- function(fit, label) {
  ll_test <- extract_log_lik(fit, "log_lik_test")
  lppd_test <- sum(log(colMeans(exp(ll_test))))
  
  yhat_test <- rstan::extract(fit, "yhat_test")$yhat_test
  yhat_test_median <- apply(yhat_test, 2, median)
  
  bind_cols(tibble(model = label),
            stan_diagnostics(fit),
            tibble(lppd_test, yhat_test = yhat_test_median))
}



### Bayesian sampling distributions ############################################
standardize_U <- function(U, j) {
  u <- U[,j]
  if (sum(u) > 1) { u <- u / sum(u) }
  return(c(u, 1 - sum(u))) # add mortality
}


posterior_array <- function(p,
                            names,
                            v1 = seq_len(dim(p)[1]),
                            v2 = seq_len(dim(p)[2]),
                            v3 = seq_len(dim(p)[3])) {
  
  out <- expand.grid(v1 = v1, v2 = v2, v3 = v3)
  out$p <- c(p)
  names(out) <- names
  return(as_tibble(out))
}


posterior_mat <- function(p,
                          names,
                          v1 = seq_len(dim(p)[1]),
                          v2 = seq_len(dim(p)[2])) {
  
  out <- expand.grid(v1 = v1, v2 = v2)
  out$p <- c(p)
  names(out) <- names
  return(as_tibble(out))
}


tr_bayes_U <- function(j, matU, N, poolU) {
  
  posU <- mean(matU) > 0  # possible tr
  mat_dim <- nrow(posU)   # number of stage classes
  
  # possible transitions for column j
  posU_col <- posU[,j]
  posU_col_id <- which(posU_col)
  
  # transition rates for relevent columns, list format
  colU <- lapply(matU, standardize_U, j = j)
  
  # column-specific sample sizes
  colN <- vapply(N, function(x) x[j], numeric(1))
  
  # for which years are transitions based on pooled data
  pooled_u <- sapply(poolU, function(x) x[j])
  
  # total number of years in dataset
  n_year_tot <- length(matU)
  
  # check NA
  check_NA_U <- all(is.na(colN))
  
  if(check_NA_U) {
    
    real_U <- sapply(matU, function(x) x[posU_col,j])
    if (!is.matrix(real_U)) real_U <- t(matrix(real_U))
    
    fit_u_raw <- real_U %>% 
      as_tibble(.name_repair = "minimal") %>% 
      setNames(seq_len(n_year_tot)) %>% 
      mutate(col = as.integer(j), row = posU_col_id) %>% 
      gather(year, tr_U, -col, -row) %>% 
      mutate(year = as.integer(year))
    
    fit_u_tidy <- tibble(rep = 1:1000, fit = list(fit_u_raw)) %>% 
      unnest() %>% 
      dplyr::select(rep, year, col, row, tr_U) %>% 
      arrange(year, rep, col, row)
    
  } else {
    
    # were transition rates estimated in given year, or based on pooled/averaged data
    year_est_u <- which(!pooled_u)
    year_pool_u <- which(pooled_u)
    
    n_year_est_u <- length(year_est_u)
    n_year_pool_u <- length(year_pool_u)
    
    # transition integer counts for all years
    colU_counts <- mapply(function(x, y) round(x * y), colU, colN)
    
    # transition integer counts for possible transitions, and directly estimated years
    colU_counts_use <- colU_counts[c(posU_col, TRUE), !pooled_u] # extra TRUE for mortality
    
    # transform count totals to individual observations
    obs_u_l <- apply(colU_counts_use, 2, function(x) rep(1:length(x), x))
    
    # create year index for each observation
    year_l <- lapply(1:length(obs_u_l), function(i) rep(i, length(obs_u_l[[i]])))
    
    # arrange data for stan
    d <- list(N = length(unlist(obs_u_l)), K = nrow(colU_counts_use),
              J = ncol(colU_counts_use), J_unobs = max(n_year_pool_u, 1),
              group = unlist(year_l), y = unlist(obs_u_l))
    
    # fit stan model
    stan_fit_u <- stanfn(stan_multinom_hier, data = d)
    
    # extract posterior samples for theta
    fit_theta <- rstan_extract(stan_fit_u, "theta")
    fit_theta_new <- rstan_extract(stan_fit_u, "theta_new")
    
    # arrange posterior samples
    df_theta <- posterior_array(fit_theta,
                           names = c("rep", "row", "year", "tr_U"),
                           v3 = year_est_u)
    
    if (n_year_pool_u > 0) {
      df_theta_new <- posterior_array(fit_theta_new,
                                       names = c("rep", "row", "year", "tr_U"),
                                       v3 = year_pool_u)
    } else {
      df_theta_new <- NULL
    }
    
    fit_u_tidy <- bind_rows(df_theta, df_theta_new) %>% 
      mutate(col = j) %>% 
      select(rep, year, col, row, tr_U) %>% 
      arrange(year, rep, col, row)
  }
  return(fit_u_tidy)
}


tr_bayes_col_F <- function(k, pooled, colN, colF_counts_use) {
  
  # which years pooled vs estimated
  year_est <- which(!pooled)
  year_pool <- which(pooled)
  
  n_year_est <- length(year_est)
  n_year_pool <- length(year_pool)
  
  # arrange data for stan
  d <- list(J = length(colF_counts_use[k,]), J_unobs = max(n_year_pool, 1),
            group = year_est, y = colF_counts_use[k,],
            offset = log(colN[!pooled]))
  
  # fit stan model
  stan_fit_f <- stanfn(stan_poisson_hier, data = d)
  
  # extract posterior samples for alpha
  fit_alpha <- rstan_extract(stan_fit_f, "alpha")
  fit_alpha_new <- rstan_extract(stan_fit_f, "alpha_new")
  
  # arrange posterior samples
  df_alpha <- posterior_mat(fit_alpha,
                            names = c("rep", "year", "tr_F"),
                            v2 = year_est)
  
  if (n_year_pool > 0) {
    df_alpha_new <- posterior_mat(fit_alpha_new,
                                  names = c("rep", "year", "tr_F"),
                                  v2 = year_pool)
  } else {
    df_alpha_new <- NULL
  }
  
  fit_f_tidy <- bind_rows(df_alpha, df_alpha_new) %>% 
    mutate(row = k) %>% 
    mutate(tr_F = exp(tr_F))
  
  return(fit_f_tidy)
}



tr_bayes_F <- function(j, matF, N, poolF, fecund) {
  
  posF <- mean(matF) > 0  # possible tr
  mat_dim <- nrow(posF)   # number of stage classes
  n_year_tot <- length(matF)  # total number of years in dataset
  
  # possible transitions for column j
  posF_col <- posF[,j]
  posF_col_id <- which(posF_col)
  
  # transition rates for relevent columns, list format
  colF <- lapply(matF, function(x) x[,j])
  
  # column-specific sample sizes
  colN <- sapply(N, function(x) x[j])
  
  if(!any(posF_col) | !fecund) {
    
    fit_f_tidy <- expand.grid(rep = 1:1000,
                              year = 1:n_year_tot,
                              col = j,
                              row = 1:mat_dim) %>% 
      as_tibble() %>% 
      mutate(tr_F = NA_real_) %>% 
      arrange(year, rep, col, row)
    
  } else {
    
    # for which years are transitions based on pooled data
    pooled_f <- sapply(poolF, function(x) x[j])

    # integer counts of new recruits
    colF_counts <- mapply(function(x, y) round(x * y), colF, colN)
    
    # integer counts for possible transitions, and directly estimated years
    colF_counts_use <- colF_counts[, !pooled_f]
    
    # fit stan models
    fit_f_tidy <- bind_rows(lapply(posF_col_id,
                                   tr_bayes_col_F,
                                   pooled = pooled_f,
                                   colN = colN,
                                   colF_counts_use = colF_counts_use)) %>% 
      mutate(col = j) %>% 
      dplyr::select(rep, year, col, row, tr_F) %>% 
      arrange(year, rep, col, row)
  }
  return(fit_f_tidy)
}


tr_bayes_wrap <- function(j, matU, matF, N, poolU, poolF, fecund = TRUE) {
  
  tr_U <- tr_bayes_U(j, matU, N, poolU)
  tr_F <- tr_bayes_F(j, matF, N, poolF, fecund)
  full_join(tr_U, tr_F, by = c("rep", "year", "col", "row"))
}



### MPM manipulations ##########################################################
MakeMat <- function(tr, mat_dim) {
  return(matrix(tr, nrow = mat_dim, ncol = mat_dim))
}


GetReproN <- function(N, posF) {
  which_repro <- apply(posF, 2, any)
  return(sum(N[which_repro]))
}


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


scale_U <- function(matU) {
  out <- apply(matU, 2, function(x) if (any(sum(x) > 1)) {x / sum(x)} else {x})
  dimnames(out) <- dimnames(matU)
  return(out)
}


mat_mean2 <- function(l, na.rm = TRUE, replace_na = TRUE) {
  m <- popbio::mean.list(l, na.rm = na.rm)
  if (replace_na) m[is.na(m)] <- 0
  return(m)
}



### MPM and age-from-stage analyses ############################################
repro_prop_start <- function(matU, start, repro_stages) {
  
  if (sum(repro_stages) == 1) {
    n <- as.numeric(repro_stages)
  } else {
    primeU <- matU
    primeU[,repro_stages] <- 0
    N <- try(solve(diag(nrow(primeU)) - primeU))
    if (class(N) == "try-error") {
      n <- NA
    } else {
      n <- rep(0, nrow(matU))
      n[repro_stages] <- N[repro_stages,start] / sum(N[repro_stages,start])
    }
  }
  return(n)
}


lx_submax <- function(lx, tmax, strip_zero = TRUE) {
  upp <- min(tmax, length(lx))
  lx <- lx[1L:upp] / lx[1L]
  if (strip_zero) lx <- lx[lx > 0]
  return(lx)
}


lx_from_mature <- function(matU, n1, crit = 0.00001, nmax = 1e4) {
  lx_vec <- numeric(1)
  lx <- 1.0
  n <- n1
  t <- 0L
  
  while (lx > crit & t < nmax) {
    n <- matU %*% n
    lx <- sum(n)
    t <- t + 1L
    lx_vec[t] <- lx
  }
  
  return(c(1, lx_vec[lx_vec > 0]))
}


qsd <- function(matU, n1, conv = 0.01, nmax = 1e4L) {

  start <- which(n1 > 0)
  
  if (!isErgodic(matU)) {
    
    nonzero <- rep(FALSE, nrow(matU))
    nonzero[start] <- TRUE
    
    n <- n1
    t <- 1L
    
    while (!all(nonzero) & t < (nrow(matU) * 3)) {
      n <- matU %*% n
      nonzero[n > 0] <- TRUE
      t <- t + 1L
    }
    
    matU <- as.matrix(matU[nonzero,nonzero])
    n1 <- n1[nonzero]
    start <- which(which(nonzero) %in% start)
  }
  
  w <- stable.stage(matU)
  n <- n1
  dist <- conv + 1
  t <- 0L
  
  while (!is.na(dist) & dist > conv & t < nmax) {
    dist <- 0.5 * (sum(abs(n - w)))
    n <- matU %*% n
    n <- n / sum(n)
    t <- t + 1L
  }
  
  return(ifelse(is.na(dist) | dist > conv, NA_integer_, t)) 
}


shape_surv2 <- function(lx, q) {
  upp <- min(q, length(lx))
  ifelse(q < 4, NA_real_, shape_surv(lx[1:upp]))
}


# stages_reached <- function(matU, rep_stages) {
#   n <- as.numeric(rep_stages)
#   nonzero <- rep_stages
#   t <- 1L
#   
#   while (!all(nonzero) & t < (nrow(matU) * 3)) {
#     n <- matU %*% n
#     nonzero[n > 0] <- TRUE
#     t <- t + 1L
#   }
# 
#   return(nonzero)
# }

