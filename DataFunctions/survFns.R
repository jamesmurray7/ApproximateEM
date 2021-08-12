#' #######
#' Functions to obtain data objects from the survival sub-model,
#' given a coxph fit, data and initial conditions for the baseline hazard (optional).
#' ---
#' These are:
#' GLOBAL:::
#' l0: hazard associated with failure times
#' ft: Failure times (non-censored)
#' SUBJECT-SPECIFIC:::
#' Fu: design matrix on all failure times survived by subject i
#' Fi: design matrix of subject i's failure time (1, Ti).
#' l0u: hazard UP TO subject i's failure time
#' l0i: hazard associated with subject i's failure time
#' ---
#' #######

surv.mod <- function(cph, data, l0.init = NULL){
  uids <- unique(data$id)
  if("coxph.null"%in%class(cph)) message("Null model")
  # Survfit
  sf <- summary(survfit(cph))
  
  # initialise empty stores
  Fi <- matrix(NA, nr = length(uids), nc = 2)
  Di <- l0i <- c()
  Fu <- l0u <- surv.times <- list()
  
  if(is.null(l0.init)){
    l0 <- diff(c(0, sf$cumhaz))
  }else{
    l0 <- l0.init
  }
  
  # loop
  for(i in uids){
    i.dat <- subset(data, id == i)
    Di[i] <- unique(i.dat$status)
    surv.times[[i]] <- which(sf$time <= unique(i.dat$survtime))
    Fi[i,] <- c(1, unique(i.dat$survtime))
    Fu[[i]] <- cbind(1, sf$time[which(sf$time <= unique(i.dat$survtime))])
    l0u[[i]] <- l0[which(sf$time <= unique(i.dat$survtime))]
    if(Di[i] == 1) l0i[i] <- l0[which(sf$time == unique(i.dat$survtime))] else l0i[i] <- 0
    # Check if censored before first failure time
    if(Di[i] == 0 & unique(i.dat$survtime) <= min(sf$time)){ l0u[[i]] <- 0; Fu[[i]] <- cbind(0, 0) }
  }
  
  nev <- c(); surv.ids <- list()
  p <- 1
  for(i in sf$time){
    nev[p] <- length(unique(data[which(data$survtime == i),]$id))
    surv.ids[[p]] <- unique(data[which(data$survtime >= i),]$id)
    p <- p+1
  }
  
  # output
  return(list(
    ft = sf$time,
    l0 = l0,
    nev = nev,
    surv.ids = surv.ids,
    surv.times = surv.times,
    l0i = l0i,
    Di = Di,
    Fi = Fi,
    Fu = Fu,
    l0u = l0u
  ))
}