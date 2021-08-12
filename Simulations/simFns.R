#' ###
#' Functions to help with arranging simulated data
#' NB: This data must come from joineRML::simData (or similar user-written fn).
#' ###

# CastDatas (K = {3, 5, 10}) ----------------------------------------------
castData3 <- function(x){
  dat <- dplyr::left_join(
    dplyr::select(x$longdat, id, cont=ctsxl, bin=binxl, time, Y.1, Y.2, Y.3),
    dplyr::select(x$survdat, id, survtime, status = cens), "id"
  ) %>% 
    mutate(status = ifelse(survtime == max(survtime), 0, status))
  
  survdat <- dplyr::distinct(dat, id, survtime, status, cont, bin)
  ph <- coxph(Surv(survtime, status)~cont + bin,survdat)
  list(dat=dat, survdat=survdat, ph=ph, 
       mst = median(survdat$survtime), prop=sum(survdat$status)/nrow(survdat))
}

castData5 <- function(x){
  dat <- dplyr::left_join(
    dplyr::select(x$longdat, id, cont=ctsxl, bin=binxl, time, Y.1:Y.5),
    dplyr::select(x$survdat, id, survtime, status = cens), "id"
  ) %>% 
    mutate(status = ifelse(survtime == max(survtime), 0, status))
  
  survdat <- dplyr::distinct(dat, id, survtime, status, cont, bin)
  ph <- coxph(Surv(survtime, status)~cont + bin,survdat)
  list(dat=dat, survdat=survdat, ph=ph, 
       mst = median(survdat$survtime), prop=sum(survdat$status)/nrow(survdat))
}

castData10 <- function(x){
  dat <- dplyr::left_join(
    dplyr::select(x$longdat, id, cont=ctsxl, bin=binxl, time, Y.1:Y.10),
    dplyr::select(x$survdat, id, survtime, status = cens), "id"
  ) %>% 
    mutate(status = ifelse(survtime == max(survtime), 0, status))
  
  survdat <- dplyr::distinct(dat, id, survtime, status, cont, bin)
  ph <- coxph(Surv(survtime, status)~cont + bin,survdat)
  list(dat=dat, survdat=survdat, ph=ph, 
       mst = median(survdat$survtime), prop=sum(survdat$status)/nrow(survdat))
}
