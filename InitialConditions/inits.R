#' ####
#' inits.R 
#' Initial conditions for longitudinal and survival sub-models
#' ---
#' Survival part STOPs here.
#' Longitudinal then undergoes MVLME fit
#' ####

# var.e, beta, D inits using lme4 -----------------------------------------
Longit.inits <- function(K, data){
  lfitK <- list()
  for(k in 1:K){
    lfitK[[k]] <- lmer(as.formula(paste0('Y.', k, '~ time + cont + bin + (1 + time|id)')), data = data)
  }
  # The three items
  var.e <- do.call(c, lapply(lfitK, sigma))^2
  beta <- do.call(c, lapply(lfitK, fixef))
  D <- as.matrix(Matrix::bdiag(
    lapply(lfitK, function(x) matrix(VarCorr(x)$id, dim(VarCorr(x)$id)))
  ))
  list(var.e.init = var.e,
       beta.init = beta,
       D.init = D,
       long.fits = lfitK)
}

# Populate REs matrix -----------------------------------------------------
Ranefs <- function(longfits){
  fits <- longfits$long.fits
  K <- length(fits)
  
  # Just need ONE instance of ID x time
  idtime.list <- list()
  for(i in unique(fits[[1]]@frame$id)) idtime.list[[i]] <- cbind(1, fits[[1]]@frame[fits[[1]]@frame$id == i, 'time'])
  
  # The random effects
  ranefK <- list()
  for(k in 1:K){
    ranefK[[k]] <- as.matrix(ranef(fits[[k]])$id)
    colnames(ranefK[[k]]) <- paste0(c("intercept_", "slope_"), k)
  }

  REs <- as.data.frame(do.call(cbind, ranefK))
  REs$id <- 1:nrow(REs)
  REs
}

# Survival Inits ----------------------------------------------------------

# Getting data into time1/time2 format...
ToStartStop <- function(data){
  this.subj <- list()
  uids <- unique(data$id)
  
  for(i in uids){
    i.dat <- data[data$id == i, c('id', 'time', 'survtime')]
    this.subj[[i]] <- cbind(
      id = i.dat$id,
      time1 = i.dat$time,
      time2 = c(i.dat$time[-1], unique(i.dat$survtime))
    )
  }
  as.data.frame(do.call(rbind, this.subj))
}

# Using this StartStop and getting timevarying coxph...
TimeVarCox <- function(data, REs, fixef.surv = c('cont', 'bin'),
                       survtime = 'survtime', status = 'status'){
  # Prepare data
  ss <- ToStartStop(data)
  ss2 <- dplyr::left_join(ss, REs, 'id')
  ss2 <- dplyr::distinct(dplyr::left_join(ss2, data[, c('id', fixef.surv, survtime, status)], 'id'))
  
  # Create \gamma_k variables
  K <- ncol(REs) %/% 2
  gammaK <- matrix(NA, nrow(ss2), K)
  colnames(gammaK) <- paste0('gamma_', 1:K)
  
  for(k in 1:K){
    ssK <- ss2[, c('time1', paste0(c('intercept_', 'slope_'), k))]
    gammaK[, k] <- ssK[, 2] + ssK[, 1] * ssK[, 3]
  }
  
  # And join on ...
  ss3 <- cbind(ss2, gammaK)
  # Update this to deal with ties too?
  ss3$status2 <- ifelse(ss3$survtime == ss3$time2, ss3$status, 0)
  
  # Time Varying coxph
  # Formula
  timevar.formula <- as.formula(
    paste0('Surv(time1, time2, status2) ~ ', paste0(fixef.surv, collapse = ' + '), ' + ', paste0('gamma_', 1:K, collapse = ' + '))
  )
  ph <- coxph(timevar.formula, data = ss3)
  
  list(inits = coef(ph), l0.init = coxph.detail(ph)$haz, ph = ph)
}
