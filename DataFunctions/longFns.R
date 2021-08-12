#' ####
#' Functions to obtain data matrices from the longitudinal sub-model
#' ---
#' Xi: Design matrix for the fixed effects - (intercept), time, cont, bin;
#' Ki: Design vector for the fixed effects in the survival model - cont, bin;
#' Yi: Outcome vector - Y.1, ..., Y.K;
#' Zi: Design matrix for the random effects - (intercept) and time;
#' mi: The number of responses for each subject for each longitudinal covariate;
#' split(X/Y/Z)ks: Splits Xi/Yi/Zi into the k = 1,...,K covariates.
#' ####

# Xi ----------------------------------------------------------------------

getXi <- function(data, K = 5){
  uids <- unique(data$id)
  xx <- list()
  for(i in uids){
    i.dat <- subset(data, id == i)
    Xi <- list()
    for(j in 1:K){
      Xi[[j]] <- cbind(1, i.dat$time, i.dat$cont, i.dat$bin)
    }
    xx[[i]] <- as.matrix(Matrix::bdiag(Xi))
    colnames(xx[[i]]) <- NULL
  }
  xx
}

# Ki ----------------------------------------------------------------------

getKi <- function(data, ph){ 
  data <- data[order(data$survtime),]
  uids <- unique(data$id)
  mm <- model.matrix(ph)
  kk <- list()
  for(i in uids){
    kk[[i]] <- matrix(mm[i, ], nr = 1)
  }
  kk
}

# Yi ----------------------------------------------------------------------

getYi <- function(data, K = 5){
  uids <- unique(data$id)
  yy <- list()
  for(i in uids){
    Yi <- c()
    for(j in 1:K) Yi <- c(Yi, data[data$id == i, paste0("Y.", j)])
    yy[[i]] <- matrix(Yi, nc = 1)
  }
  yy
}

# Zi ----------------------------------------------------------------------

getZi <- function(data, K = 5){
  uids <- unique(data$id)
  zz <- list()
  for(i in uids){
    Zi <- list()
    Zij <- cbind(1, data[data$id == i, "time"])
    for(j in 1:K) Zi[[j]] <- Zij
    zz[[i]] <- as.matrix(Matrix::bdiag(Zi))
  }
  zz
}

# mi ----------------------------------------------------------------------

getmi <- function(data, K = 5){
  uids <- unique(data$id)
  mm <- list()
  for(i in uids){
    mi <- c()
    for(j in 1:K) mi <- c(mi, length(data[data$id == i, paste0("Y.", j)]))
    mm[[i]] <- mi
  }
  mm
}

# Split Y, X, Z into Yk, Xk, Zk -------------------------------------------

splitYks <- function(Y, mi, K = 5){
  uids <- length(Y) # a list
  kk <- list()
  for(i in 1:uids){
    ll <- Y[[i]]; mm <- c(0, cumsum(mi[[i]]))
    Yik <- list()
    for(j in 2:(K+1)) Yik[[j-1]] <- ll[(mm[j-1] + 1):mm[j]]
    kk[[i]] <- do.call(cbind, Yik)
  }
  kk
}

splitXks <- function(data, K = 5){
  uids <- unique(data$id)
  kk <- list()
  for(i in uids){
    Xik <- list() # predicated on balanced design.
    Xi <- cbind(1, data[data$id == i, c("time", "cont", "bin")])
    for(j in 1:K) Xik[[j]] <- as.matrix(Xi)
    kk[[i]] <- Xik
  }
  kk
}

splitZks <- function(data, K = 5){
  uids <- unique(data$id)
  kk <- list()
  for(i in uids){
    Zik <- list()
    Zi <- cbind(1, data[data$id == i, "time"])
    for(j in 1:K) Zik[[j]] <- as.matrix(Zi)
    kk[[i]] <- Zik
  }
  kk
}
