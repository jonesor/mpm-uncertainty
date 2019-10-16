

### libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(popbio)
library(Rage)
library(tidyverse)
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


