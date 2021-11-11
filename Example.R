#' #####
#' Example.R        // Assumes working directory is ~/.../ApproximateEM/
#' --------
#' Section 1 of this program illustrates the simulation of (trivariate) data under the joint modelling framework
#' using the simData function from R package joineRML. And subsequent joint model fit to the data using the 
#' approximate EM algorithm outlined in Murray & Philipson's 
#' "A Fast Approximate EM Algorithm for Joint Models of Survival and Multivariate Longitudinal Data".
#' --------
#' Section two of this program simulates a set of N data sets and iteratively applies the approximate
#' EM Algorithm to each of the N data sets, before returning a list of model fits along with an example of inspection.
#' --------
#' The simulations performed here use th same true parameter estimates as section 3.1 in the paper,
#' which in turn are the same as used in earlier work by Philipson et al. (2020) (also referenced in the main paper).
#' #####
rm(list=ls())
source('./EM.R')

# Defining a helper function to take the output of joineRML::simData and make it usable by em().
castData <- function(x, nK){
  dat <- dplyr::left_join(
    dplyr::select(x$longdat, id, cont=ctsxl, bin=binxl, time, paste0('Y.', 1:nK)),
    dplyr::select(x$survdat, id, survtime, status = cens), "id"
  ) %>% 
    mutate(status = ifelse(survtime == max(survtime), 0, status))
  
  survdat <- dplyr::distinct(dat, id, survtime, status, cont, bin)
  ph <- coxph(Surv(survtime, status)~cont + bin,survdat)
  list(dat=dat, survdat=survdat, ph=ph)
}

# Section 1 ---------------------------------------------------------------
# Define true parameters
beta <- rbind(c(0, 1, 1, 1),
              c(0, -1, 0, 0.5),
              c(0, 0, 0.5, -0.5))

D <- diag(6)
D[1, 1] <- D[3, 3] <- D[5, 5] <- 0.5^2
D[2, 2] <- D[4, 4] <- D[6, 6] <- 0.2^2
D[1, 3] <- D[3, 1] <- -0.5 * 0.5 * 0.5
D[1, 5] <- D[5, 1] <- 0.5^3
D[3, 5] <- D[5, 3] <- -0.5*(0.5^2)

var.e <- c(0.25, 0.25, 0.25)
eta <- c(0, 1)
gamma <- c(-0.5, 1, 0.5)

# Simulate some data and transform it using the above castData function.
sim <- joineRML::simData(n = 250, ntms = 10, beta = beta, gamma.x = eta, gamma.y = gamma, sigma2 = var.e,
                         D = D, theta0 = -4, theta1 = 0.15)

dat <- castData(sim, 3)

# Fit using EM(data, ph, gh.nodes, nK)
fit <- em(dat$dat, dat$ph, gh.nodes = 3, nK = 3)

# And inspect coefficients and their SEs (random ones chosen)
source('Simulations/simFns.R')
extract.estimates(fit, 3)
fit$SEs

# Section 2 ---------------------------------------------------------------
rm(beta, D, var.e, gamma, eta, sim, dat, fit)

# Define number of simulated data sets. 10 obviously far too small and is for demonstration purposes only.
N <- 10 
# Define true parameters
beta <- rbind(c(0, 1, 1, 1),
              c(0, -1, 0, 0.5),
              c(0, 0, 0.5, -0.5))

D <- diag(6)
D[1, 1] <- D[3, 3] <- D[5, 5] <- 0.5^2
D[2, 2] <- D[4, 4] <- D[6, 6] <- 0.2^2
D[1, 3] <- D[3, 1] <- -0.5 * 0.5 * 0.5
D[1, 5] <- D[5, 1] <- 0.5^3
D[3, 5] <- D[5, 3] <- -0.5*(0.5^2)

var.e <- c(0.25, 0.25, 0.25)
eta <- c(0, 1)
gamma <- c(-0.5, 1, 0.5)

# Generate N samples
sim <- replicate(N,
                 joineRML::simData(n = 250, ntms = 10, beta = beta, gamma.x = eta, gamma.y = gamma, sigma2 = var.e,
                                   D = D, theta0 = -4, theta1 = 0.15),
                 simplify = F)


dat <- lapply(sim, castData, 3)

# Fit N models with a handy progress bar(!)
fits <- list()
pb <- utils::txtProgressBar(max = N, style = 3)
for(z in 1:N){
  fits[[z]] <- suppressMessages(
    tryCatch(em(dat[[z]]$dat, dat[[z]]$ph, gh.nodes = 3, nK = 3), 
             error = function(e) NULL)
  )
  utils::setTxtProgressBar(pb, z)
}

# Extract coefficients of interest (say, \gamma and \eta)
gamma.eta <- lapply(fits, function(x) c(x$coeffs$gamma, x$coeffs$eta))
apply(do.call(rbind, gamma.eta), 2, mean)

# And see how the approximate SE compares to empirical SD
gamma.eta.SE <- lapply(fits, function(x){
  ses <- x$SEs
  ses[grepl('gamma|^eta', names(ses))]
})
apply(do.call(rbind, gamma.eta), 2, sd)
apply(do.call(rbind, gamma.eta.SE), 2, mean)

# or indeed extract all
do.call(rbind, lapply(fits, extract.estimates, 3))

# Get some indication of elapsed time
elapsed <- do.call(c, lapply(fits, function(x) x$EMtime + x$postprocess.time))
quantile(elapsed, c(.25, .5, .75))
