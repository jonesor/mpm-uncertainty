

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
string_to_mat <- function(A) {
  A <- gsub(pattern = "\\[|\\]|\\;", "", A)
  A <- strsplit(x = A, split = " ")[[1]]
  mat <- matrix(as.numeric(A), nrow = sqrt(length(A)), byrow = TRUE)
  return(mat)
}


nx_to_vec <- function(x) {
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

rdata_load2 <- function(path) {
  env <- new.env()
  x <- load(path, env)[1]
  out <- env[[x]]
  out$Altitude <- NULL
  out$MatrixStartYear <- NULL
  out$MatrixEndYear <- NULL
  return(out)
}



### Confirm that stage-specific sample sizes match transition rates ############
check_freqs_mat <- function(matU, N, prec = 0.001) {
  dim <- nrow(matU)
  out <- character(dim)
  
  for (i in 1:dim) {
    out[i] <- check_freqs(N[i], matU[,i], prec = prec) 
  }
  
  return(data.frame(N, x = out))
}

check_freqs <- function(n, x, prec = 0.001) {
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


sim_stage_U <- function(x, vital_ind, n) {
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


sim_stage_F <- function(x, vital_ind, n) {
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


sim_U <- function(matU, posU, N) {
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
      simU[,i] <- sim_stage_U(matU[,i], vital_ind, N[i])
    }
  }
  return(simU)
}


sim_F <- function(matF, posF, N) {
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
      simF[,i] <- sim_stage_F(matF[,i], vital_ind, N[i])
    }
  }
  return(simF)
}


sim_U_wrapper <- function(matU, posU = matU > 0, N, nsim) {
  return(replicate(nsim, sim_U(matU, posU, N), simplify = FALSE))
}


sim_F_wrapper <- function(matF, posF = matF > 0, N, nsim) {
  return(replicate(nsim, sim_F(matF, posF, N), simplify = FALSE))
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



### MPM manipulations ##########################################################
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
lx_submax <- function(lx, tmax, strip_zero = TRUE) {
  upp <- min(tmax, length(lx))
  lx <- lx[1L:upp] / lx[1L]
  if (strip_zero) lx <- lx[lx > 0]
  return(lx)
}


shape_surv2 <- function(lx, q) {
  upp <- min(q, length(lx))
  ifelse(q < 4, NA_real_, shape_surv(lx[1:upp]))
}


sum2 <- function(x) {
  ifelse(all(is.na(x)), NA_real_, sum(x, na.rm = TRUE))
}


pool_counts <- function(nl) {
  X <- do.call(rbind, nl)
  return(apply(X, 2, sum2))
}


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
  posU[exclude,] <- FALSE
  posU[,exclude] <- FALSE
  
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



