#' ###
#' Functions to help with arranging simulated data and extracting from fitted objects.
#' NB: This data must come from joineRML::simData (or similar user-written fn).
#' ###

# CastDatas (K = {3, 5, 10}) ----------------------------------------------
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

# Extracting results ------------------------------------------------------
extract.estimates <- function(x, nK){ # x a fitted list.
  # Covariance matrix D
  D <- x$coeffs$D; vD <- vech(D)
  names(vD) <- paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste0, collapse = ','),']')
  # beta
  beta <- c(x$coeffs$beta); names(beta) <- names(x$SEs[grepl('^beta', names(x$SEs))])
  # sigma^2
  var.e <- c(x$coeffs$var.e); names(var.e) <- paste0('var.e_', 1:nK)
  # gamma, eta
  gamma <- c(x$coeffs$gamma)
  eta <- c(x$coeffs$eta); names(eta) <- paste0('eta_', names(eta))
  return(c(vD, beta, var.e, gamma, eta))
}

extract.time <- function(x) x$EMtime + x$postprocess.time
