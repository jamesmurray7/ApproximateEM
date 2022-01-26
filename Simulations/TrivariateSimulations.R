#' ####
#' TrivariateSimulations.R
#' Simulating trivariate data used in Simulation I
#' ----
#' Script loops through the profiles used in Section 3.1 in the paper, 
#' although this can be escapes and specific ones chosen via
#' joineRML::simData(.) or below simn with bespoke choices.
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
gamma.x <- c(-0.1, 0.3)
gamma.y <- c(-0.5, 0.75, 0.5)

simn <- function(n, ntms, theta){
  x <- replicate(100,
                 joineRML::simData(n = n, ntms = ntms, beta = beta, gamma.x = gamma.x, # Medium profile
                                   gamma.y = gamma.y, sigma2 = sigma2, D = D, theta0=theta[1], theta1=theta[2]), #theta0 = -4, theta1 = 0.15),
                 simplify = F)
  print(mean(do.call(c, lapply(x, function(xx) sum(xx$survdat$cens)/n))))
  x
}

source("./Simulations/simFns.R")
save.loc <- paste0("./Simulations") # By default will simply 'dump' into the Simulations folder.

for(i in c(250, 500, 1e3)){
  for(j in c("Short", "Medium", "Long")){
    if(j == "Short"){ntms <- 6; theta <- c(-2.5, 0.15)}
    if(j == "Medium"){ntms <- 10; theta <- c(-3.5, 0.15)}
    if(j == "Long"){ntms <- 15; theta <- c(-5, 0.15)}
    d <- simn(i, ntms, theta)
    dat <- lapply(d, castData, 3)
    file.name <- paste0(save.loc, "Trivariate-n", i, "-", j, ".RData")
    save(dat, file = file.name)
    message(file.name)
  }
}

# Fit using approximate EM algorithm --------------------------------------
# Load e.g. dat250 and fit EM
load("./Simulations/Trivariate-n250-Medium.RData")
onefit <- em(dat[[1]]$dat, dat[[1]]$ph)

# Or define a function and lapply
emfit <- function(x, ...){
  # print(x$aa)
  res <- tryCatch(em(x$dat, x$ph, ...), error = function(e) NULL)
  res
}

data.dir <- './Simulations/'    # Change to match above.
data.files <- as.list(dir(data.dir, pattern = '^Tr'))

fitter <- function(x){
  load(paste0(data.dir, x))
  message(paste0(data.dir, x))
  pb <- utils::txtProgressBar(max = 100, style = 3)
  fits <- list()
  for(i in 1:100){
    fits[[i]] <- suppressMessages(emfit(dat[[i]]))
    utils::setTxtProgressBar(pb, i)
  }
  save(fits, file = paste0(data.dir, 'fits-', x)) # Saves back into whatever folder chosen above.
  message(paste0(data.dir, 'fits-', x,'\n'))
  fits      
}

all.fits <- lapply(data.files, fitter)

