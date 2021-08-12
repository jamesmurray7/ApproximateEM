#' #####
#' K5-sims.R
#' Simulating data for K=5 case.
#' ----
#' The parameter values are set using ADNI (See Appendix for further details).
#' for n=500 and approximately 50% failure rate.
#' #####

# Longitudinal parameters -------------------------------------------------
# ADAS13
beta1 <- c(-0.1, 0.2, 0.1, 0.05)
var1 <- 0.2
D1 <- matrix(c(0.6, 0.1, 0.1, 0.05), 2, 2)

# FAQ
beta2 <- c(-0.1, 0.25, 0.005, -0.1)
var2 <- 0.2
D2 <- matrix(c(0.6, 0.05, 0.05, 0.1), 2, 2)

# MidTemp
beta3 <- c(0.25, -0.1, -0.2, -0.8)
var3 <- 0.05
D3 <- matrix(c(0.8, 0.025, 0.025, 0.01), 2, 2)

# RAVLT 
beta4 <- c(-0.2, -0.1, -0.1, 0.4)
var4 <- 0.2
D4 <- matrix(c(0.7, 0.05, 0.05, 0.02), 2, 2)

# Hippocampus 
beta5 <- c(0.2, -0.15, -0.4, -0.5)
var5 <- 0.01
D5 <- matrix(c(0.7, 0.02, 0.02, 0.01), 2, 2)

# Populating D, the covariance matrix -------------------------------------
D <- as.matrix(Matrix::bdiag(D1, D2, D3, D4, D5))
D[3:4, 1:2] <- matrix(c(0.2, 0.05, 0.05, 0.02), 2, 2)
D[5:6, 1:2] <- matrix(c(-0.2, -0.05, -0.01, -0.005), 2, 2)
D[7:8, 1:2] <- matrix(c(-0.5, -0.05, -0.05, -0.01), 2, 2)
D[9:10, 1:2] <- matrix(c(-0.25, -0.05, -0.01, -0.005), 2, 2)
D[5:6, 3:4] <- matrix(c(-0.1, -0.05, -0.005, -0.005), 2, 2)
D[7:8, 3:4] <- matrix(c(-0.2, -0.05, -0.1, -0.005), 2, 2)
D[9:10, 3:4] <- matrix(c(-0.15, -0.05, -0.005, -0.004), 2, 2)
D[7:8, 5:6] <- matrix(c(0.2, 0.02, 0.02, 0.002), 2, 2)
D[9:10, 5:6] <- matrix(c(0.4, 0.02, 0.01, 0.001), 2, 2)
D[9:10, 7:8] <- matrix(c(0.2, 0.05, 0.05, 0.01), 2, 2)

# Make symmetric, ensure semi-positive definite
D[upper.tri(D)] <- t(D)[upper.tri(D)] 
any(eigen(D)$values < 0)

D.PD <- as.matrix(Matrix::nearPD(D, keepDiag = T)$mat)

# Survival parameters
eta <- c(0.05, -0.3)
gamma <- c(0.6, 0.5, -0.2, -0.25, -0.2)


# Simulate data -----------------------------------------------------------
N <- 100
x <- replicate(N, 
               joineRML::simData(n = 500, ntms = 15,
                                 beta = rbind(beta1, beta2, beta3, beta4, beta5),
                                 sigma2 = c(var1, var2, var3, var4, var5),
                                 gamma.x = eta,
                                 gamma.y = gamma,
                                 D = D.PD, theta0 = -3.15, theta1 = 0.15
               ),
               simplify = F)


# checking failure rate

propfail <- function(x) sum(x$survdat$cens)/500 * 100
mean(do.call(c, lapply(x, propfail)))  # approximately 50 %.

source("./Simulations/simFns.R")
dat <- lapply(x, castData5)
