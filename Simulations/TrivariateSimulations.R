#' ####
#' TrivariateSimulations.R
#' Simulating trivariate data used in simulation I
#' ----
#' Shown is the 'medium' profile.
#' ####

beta <- rbind(c(0, 1, 1, 1),
              c(0, -1, 0, 0.5),
              c(0, 0, 0.5, -0.5))
D <- diag(6)
D[1, 1] <- D[3, 3] <- D[5, 5] <- 0.5^2
D[2, 2] <- D[4, 4] <- D[6, 6] <- 0.2^2
D[1, 3] <- D[3, 1] <- -0.5 * 0.5 * 0.5
D[1, 5] <- D[5, 1] <- 0.5^3
D[3, 5] <- D[5, 3] <- -0.5*(0.5^2)

sigma2 <- c(0.25, 0.25, 0.25)
gamma.x <- c(0, 1)
gamma.y <- c(-0.5, 1, 0.5)

simn <- function(n){
      x <- replicate(100,
          joineRML::simData(n = n, ntms = 10, beta = beta, gamma.x = gamma.x,  # Medium profile
                            gamma.y = gamma.y, sigma2 = sigma2, D = D, theta0 = -4, theta1 = 0.15),
          simplify = F)
      print(mean(do.call(c, lapply(x, function(xx) sum(xx$survdat$cens)/n))))
      x
}

source("./Simulations/simFns.R")

for(i in c(250, 500, 1000)){
  d <- simn(i)
  dat <- lapply(d, castData3)
  for(j in 1:100) dat[[j]]$aa <- j
  save(dat, file = paste0(getwd(), "/Simulations/Trivariate", i, ".RData"))
}


# Fit using approximate EM algorithm --------------------------------------

# Load e.g. dat250 and fit EM
load("./Simulations/Trivariate250.RData")
onefit <- em(dat[[1]]$dat, dat[[1]]$ph)

# Or define a function and lapply
emfit <- function(x, ...){
  print(x$aa)
  res <- tryCatch(em(x$dat, x$ph, ...), error = function(e) NULL)
  res
}

fits <- lapply(dat, emfit)
