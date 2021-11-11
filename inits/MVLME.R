#' ###
#' MVLME.R // Multivariate linear mixed effects
#' ---
#' Uses longit initial conditions from inits.R
#' and fits MVLME to find optimal starting values for joint model
#' ###

if(!'Rcpp'%in%loadedNamespaces()) library(Rcpp)
if(!'RcppArmadillo'%in%loadedNamespaces()) library(RcppArmadillo)

# Prerequisites -----------------------------------------------------------
# sourceCpp("./mvlmecpp.cpp") # Loading Rcpp Source files...

if(!"vech"%in%ls()) vech <- function(x) x[lower.tri(x, diag = T)]

splitbbT <- function(bbT, K){
  uids <- length(bbT)
  out <- list()
  for(i in 1:uids){
    out[[i]] <- lapply(split(seq(K * 2), rep(1:K, each = 2)), function(x) bbT[[i]][x, x])
  }
  out
}

mvlme <- function(d, Y, X, Z, Yk, Xk, Zk, mi, inits.long, K,
                  tol.mvlme = 5e-3){
  # Data
  diff <- 100; iter <- 0;
  uids <- unique(d$id)
  n <- length(uids)
  XtX <- lapply(X, crossprod)
  
  # Unpack inits.long
  D <- inits.long$D.init
  var.e <- inits.long$var.e.init
  beta <- inits.long$beta.init; beta.mat <- matrix(beta, nr = K, byrow = T)
  V <- lapply(mi, function(i) {
    diag(x = rep(var.e, i), ncol = sum(i))
  })
  
  # And make parameter vector
  params <- c(vech(D), var.e, beta)
  names(params) <- c(rep("D", length(vech(D))),
                     paste0("var.e_", 1:K), paste0(rep(paste0("beta_", 1:K), each=4), 0:3))
  
  # EM Algorithm
  mvlmetime <- c()
  EMstart <- proc.time()[3]
  
  while(diff > tol.mvlme){
    Estart <- proc.time()[3]
    #' E-step -------------------
    # E[b]
    b <- Eb(Y, X, Z, V, D, beta, n)
    bmat <- lapply(b, function(x) matrix(x, nc = 2, byrow = T))
    
    # Sigmai
    Sigmai <- covb(Z, V, solve(D), n)
    
    # E[bbT]
    bbT <- EbbT(b, Sigmai, n)
    bbTk <- splitbbT(bbT, K)
    
    #' M-step -------------------
    
    # D //
    D.new <- Reduce('+', bbT) / length(uids)
    
    # \beta //
    beta.rhs <- betaRHS(X, Y, Z, b, n)
    beta.new <- solve(Reduce('+', XtX)) %*% Reduce('+', beta.rhs)
    names(beta.new) <- names(beta)
    
    # var.e
    # change to beta.new.mat?
    var.e.new <- Ee(Yk, Xk, Zk, beta.mat, bmat, bbTk, n, K)/colSums(do.call(rbind, mi))
    
    mvlmetime[iter + 1] <- proc.time()[3] - Estart
    
    # New parameter vector
    params.new <- c(vech(D.new), var.e.new, beta.new)
    names(params.new) <- names(params)
    # Take difference & report
    diffs <- abs(params.new - params)/(abs(params) + 1e-3)
    diff <- max(diffs)
    message("\nIteration ", iter + 1, " largest relative difference = ", round(diff, 5))
    message("For: ", names(params)[which(diffs == diff)])
    message("---> Old: ", params[which(diffs == diff)], ", new: ", params.new[which(diffs == diff)])
    
    # Update parameters and loop
    params <- params.new
    D <- D.new; var.e <- var.e.new; 
    beta <- beta.new; beta.mat <- matrix(c(beta), nr = K, byrow = T)
    V <- lapply(mi, function(i) {
      diag(x = rep(var.e.new, i), ncol = sum(i))
    })
    iter <- iter + 1
  }
  message("Converged after ", iter, " iterations ")
  list(
    beta = beta, D = D, var.e = var.e, b = b,
    V = V, XtX = XtX,
    elapsed.time = proc.time()[3] - EMstart, EMtime = sum(mvlmetime)
  )
}
