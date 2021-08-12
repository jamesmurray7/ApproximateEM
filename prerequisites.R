#' ####
#' Prerequisites.R
#' ----
#' Loads all functions and libraries necessary for EM.R (trivariate)
#' ####

# Libraries ---------------------------------------------------------------
library(survival)
library(lme4)
library(Rcpp)
library(RcppArmadillo)
library(dplyr)

# Cpp files ---------------------------------------------------------------
message("Loading ", paste0(getwd(), "/source/mvlme.cpp"))
sourceCpp("./source/mvlme.cpp")
message("Loading ", paste0(getwd(), "/source/gammaCalc.cpp"))
sourceCpp("./source/gammaCalc.cpp")
message("Loading ", paste0(getwd(), "/source/ll.cpp"))
sourceCpp("./source/ll.cpp")

# Source functions --------------------------------------------------------
source("./DataFunctions/longFns.R")
source("./DataFunctions/survFns.R")
source("./InitialConditions/inits.R")
source("./InitialConditions/MVLME.R")

# Other functions ---------------------------------------------------------
vech <- function(x) x[lower.tri(x, diag = T)]
tr <- function(x) sum(diag(x))
an <- as.numeric
repCols <- function(X, n = 3) X[,rep(1:2, n)] 
repVec <- function(x, n = 3) rep(x, n)

d2b.ll <- function(K, Z, D, V, l0u, Fu, g, eta, bi){
  gr <- rep(g, each = 2)
  b <- bi
  diag.g <- diag(gr)
  if(nrow(Fu) == 1){
    surv.part <- -diag.g %*% repCols(Fu) %*% (l0u * exp(K %*% eta + repCols(Fu) %*% (gr * b))) %*% repCols(Fu) %*% diag.g
  }else{
    diag.kern <- diag(l0u * an(exp(K %*% eta) %x% exp(repCols(Fu) %*% (gr * b))))
    surv.part <- -diag.g %*% crossprod(repCols(Fu), diag.kern) %*% repCols(Fu) %*% diag.g
  }
  -crossprod(Z, solve(V) %*% Z)- solve(D) + surv.part
}

message("\nPrerequisites Loaded\n")

