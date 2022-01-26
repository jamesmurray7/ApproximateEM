#' ####
#' Expectation Maximisation 
#' Assumes working directory is set to the parent directory
#'                    (i.e. ~/.../ApproximateEM)
#' ----
#' The call 'source("./prerequisites.R")' loads all functions necessary
#'       NB: Part of this includes loading several .cpp files using Rcpp
#'           and RcppArmadillo, which may require separate installation of 
#'           appropriate compiler tools depending on OS.
#'           This step can take some time in the first instance of R session.
#' 
#' The functions used to set-up data matrices require data with a column for subject ID
#' called 'id'; time variable 'time'; continuous variable 'cont' and binary 'bin'.
#' The K longitudinal outcomes should be named 'Y.1, ..., Y.K'. Currently
#' implementation is predicated upon a balanced design only.
#'  
#' The EM call itself requires input of the dataset and a coxph model fit on the 
#' 'cont' and 'bin' variables only (castData prepares data and ph in a suitable fashion from simulated data)
#' 
#' The user can control the number of Gauss-Hermite quadrature nodes used (gh.nodes); the
#' maximum number of iterations to try before exiting (max.iter); the tolerance used
#' by the algorithm (tol); the convergence criterion to use (diff.type, must be 'abs' or
#' 'abs.rel'); if parameter estimates should be stored after each iteration (collect.hist);
#'  nK denotes the number of covariates and should automate the function.
#' ####

source("./prerequisites.R")

em <- function(data, ph,
               gh.nodes = 3, collect.hist = T, max.iter = 200, 
               tol=0.01, diff.type = "abs.rel", nK = 3, post.process = T, verbose = F){
  # Set-up ------------------------------------------------------------------
  start.time <- proc.time()[3]
  q <- nK * 2
  diff <- 100; b.diff <- 100; iter <- 0
  uids <- unique(data$id); n <- length(uids)
  #' Data matrices ----
  X <- getXi(data, nK); Z <- getZi(data, nK); Y <- getYi(data, nK)
  mi <- getmi(data, nK)
  XtX <- lapply(X, crossprod)
  Zk <- splitZks(data, nK)
  Xk <- splitXks(data, nK)
  Yk <- splitYks(Y, mi, nK)
  K <- getKi(data, ph)
  
  #' Initial conditions -----
  inits.long <- suppressWarnings(Longit.inits(nK, data))
  inits.surv <- TimeVarCox(data, Ranefs(inits.long))
  
  # MVLME for optimal initial conditions given observed data for longitudinal part
  mvlme.fit <- mvlme(data, Y, X, Z, Yk, Xk, Zk, mi, inits.long, nK)
  
  #' Survival-related objects ----
  sv <- surv.mod(ph, data, inits.surv$l0.init)
  ft <- sv$ft; nev <- sv$nev
  surv.ids <- sv$surv.ids; surv.times <- sv$surv.times
  Di <- sv$Di; Deltai.list <- as.list(Di)
  l0 <- sv$l0; l0i <- sv$l0i; l0i.list <- as.list(l0i); l0u <- sv$l0u
  Fi <- sv$Fi; Fu <- sv$Fu
  # List-versions of several data objects for use of mapply()
  Fi.list <- lapply(1:nrow(Fi), function(i) Fi[i, ])
  rvFi.list <- lapply(1:nrow(Fi), function(i) do.call(c, replicate(nK, Fi[i,], simplify = F)))
  KK <- sapply(1:n, function(x){  # For updates to \eta
    x <- apply(K[[x]], 2, rep, nrow(Fu[[x]]))
    if("numeric" %in% class(x)) x <- t(as.matrix(x))
    x
  }) 
  
  # Extract initial conditions
  D <- mvlme.fit$D
  b <- lapply(mvlme.fit$b, c)
  beta <- c(mvlme.fit$beta)
  var.e <- mvlme.fit$var.e
  V <- mvlme.fit$V
  mvlme.time <- mvlme.fit$elapsed.time
  gamma <- inits.surv$inits[3:length(inits.surv$inits)]; gr <- rep(gamma, each = 2)  # 2 for intercept and slope + proportional assoc
  eta <- inits.surv$inits[1:2]
  
  # Cast to parameter vector
  vD <- vech(D)
  names(vD) <- paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste0, collapse = ','),']')
  params <- c(vD, beta, var.e, gamma, eta)
  names(params) <- c(names(vD), 
                     names(beta), names(inits.long$var.e.init), names(gamma), paste0('eta_', names(eta)))
  
  # Collect data objects and iteration "0" history
  dmats <- list(Y = Y, X = X, Z = Z, mi = mi, Yk = Yk, Xk = Xk, Zk = Zk,        # Longit.
                K = K, KK = KK, Fu = Fu, Fi = Fi.list, Deltai = Deltai.list)    # Survival
  if(collect.hist) iter.hist = data.frame(iter = iter, t(params))
  
  # Gaussian Quadrature -----------------------------------------------------
  gh <- statmod::gauss.quad.prob(gh.nodes, "normal")
  v <- gh$n; w <- gh$w
  
  # Define indices for helping later
  b.inds <- split(seq(nK * 2), rep(1:nK, each = 2))
  beta.inds <- split(seq(length(beta)), rep(1:nK, each = 4))
  
  message("Starting EM Algorithm")
  EM.time <- c()
  # EM ----------------------------------------------------------------------
  while(diff > tol & iter < max.iter){
    p1 <- proc.time()[3]
    
    #' ##################
    #' E-step ========= #
    #' ##################

    b.hat <- mapply(function(b, Y, X, Z, V, mi, K, Delta, l0i, Fi, l0u, Fu,
                             rvFi){
      ucminf::ucminf(b, ll, gradll, 
                     Y, X, Z, V, D, sum(mi), K, Delta, l0i, Fi, l0u, Fu,
                     gamma, beta, eta, gr, rvFi, nK, q,
                     control = list(xtol = 1e-3, grtol = 1e-6))$par
    },b = b, Y = Y, X = X, Z = Z, V = V, mi = mi, K = K, Delta = Deltai.list, l0i = l0i.list,
      Fi = Fi.list, l0u = l0u, Fu = Fu, rvFi = rvFi.list, SIMPLIFY = F)
    b.hat.split <- lapply(b.hat, function(y) lapply(b.inds, function(x) y[x]))
    bmat <- lapply(b.hat, matrix, nc = 2, byr = T)
    
    Sigmai <- mapply(function(b, Z, V, K, l0u, Fu){
      solve(-1 * sdll(b, Z, D, V, K, l0u, Fu, eta, gr, nK))
    }, b = b.hat, Z = Z, V = V, K = K, l0u = l0u, Fu = Fu, SIMPLIFY = F)

    S <- lapply(Sigmai, function(y) lapply(b.inds, function(x) y[x,x]))   # Split out into K constituent sub-matrices along block diagonal
    
    #' Step to update D ----
    D.newi <- mapply(function(S, b){
      S + tcrossprod(b)
    }, S = Sigmai, b = b.hat, SIMPLIFY = F)
    
    #' Steps to update longitudinal parameters ----
    #' Define \tau and \mu
    tau.long <- mapply(function(S, Z){
      sqrt(diag(tcrossprod(Z %*% S, Z)))
    }, S = Sigmai, Z = Z, SIMPLIFY = F)
    
    tau.longK <- mapply(function(S, Z){
      out <- list()
      for(k in 1:nK) out[[k]] <- sqrt(diag(tcrossprod(Z[[k]] %*% S[[k]], Z[[k]])))
      out
    }, S = S, Z = Zk, SIMPLIFY = F)
    
    mu.longK <- mapply(function(X, Z, b){
      out <- list()
      for(k in 1:nK) out[[k]] <- X[[k]] %*% beta[beta.inds[[k]]] + Z[[k]] %*% b[b.inds[[k]]]
      out
    }, X = Xk, Z = Zk, b = b.hat , SIMPLIFY = F)
    
    #' Update to \beta ----
    beta.rhs <- mapply(function(X, Y, Z, b, tau){
      rhs <- 0
      mu <- Z %*% b
      for(l in 1:gh.nodes) rhs <- rhs + w[l] * tau * v[l]
      crossprod(X, Y - mu - rhs)
    }, X = X, Y = Y, Z = Z, b = b.hat, tau = tau.long, SIMPLIFY = F)
    
    #' Update to \sigma^2_\epsilon ----
    Ee <- mapply(function(Y, mu, tau){
      temp <- matrix(NA, nr = gh.nodes, nc = nK)
      for(k in 1:nK){
        for(l in 1:gh.nodes){
          temp[l, k] <- w[l] * crossprod(Y[, k] - mu[[k]] - tau[[k]] * v[l])
        }
      }
      colSums(temp)
    }, Y = Yk, mu = mu.longK, tau = tau.longK, SIMPLIFY = F)
    
    #' #####
    #' Survival Parameters
    #' #####
    
    #' Set out Newton-Raphson items for update to (\gamma, \eta)
    Sge <- mapply(function(bmat, S, K, KK, Fu, Fi, l0u, Delta){
      Sgammaeta(c(gamma, eta), bmat, S, K, KK, Fu, Fi, l0u, Delta, w, v, 1e-4)
    }, bmat = bmat, S = S, K = K, KK = KK, Fu = Fu, Fi = Fi.list, 
    l0u = l0u, Delta = Deltai.list, SIMPLIFY = F)
    
    Hge <- mapply(function(bmat, S, K, KK, Fu, Fi, l0u, Delta){
      Hgammaeta(c(gamma, eta), bmat, S, K, KK, Fu, Fi, l0u, Delta, w, v, 1e-4)
    }, bmat = bmat, S = S, K = K, KK = KK, Fu = Fu, Fi = Fi.list, 
    l0u = l0u, Delta = Deltai.list, SIMPLIFY = F)
    
    #' ##################
    #' M-step ========= #
    #' ##################
    
    #' D -----
    D.new <- Reduce('+', D.newi)/n
    
    #' \beta ----
    beta.new <- solve(Reduce('+', XtX)) %*% Reduce('+', beta.rhs) # NB this slightly faster than Reduce('+',.) on rhs.
    #' var.e ----
    var.e.new <- colSums(do.call(rbind, Ee))/colSums(do.call(rbind, mi))      # NB this same speed as storing Ee directly as an array.
    #' The baseline hazard, \lambda ----
    lambda <- lambdaUpdate(surv.times, ft, gamma, eta, K, S,
                           b.hat.split, n, w, v, gh.nodes, nK)                # NB this not monitored for convergence.
    l0.new <- nev/rowSums(lambda)
    l0u.new <- lapply(l0u, function(x){
      ll <- length(x); l0.new[1:ll]
    })
    l0i.new <- c()
    l0i.new[which(Di == 0)] <- 0 
    l0i.new[which(Di == 1)] <- l0.new[match(Fi[which(Di==1), 2], ft)]
    
    #' (\gamma, \eta) ----
    gamma.eta.new <- c(gamma, eta) - solve(Reduce('+', Hge), rowSums(do.call(cbind, Sge)))
    
    gamma.new <- gamma.eta.new[1:nK]
    eta.new <- gamma.eta.new[(nK + 1):length(gamma.eta.new)]
    
    EM.time[iter + 1] <- proc.time()[3] - p1 # M step finishes here ---
    
    #' Update parameters and print ----
    params.new <- c(vech(D.new), beta.new, var.e.new, gamma.new, eta.new); names(params.new) <- names(params)
    if(verbose) print(sapply(params.new, round, 4))
    # Take differences (user input)
    if(diff.type == "abs"){
      diffs <- abs(params.new-params)
      b.diff <- max(abs(do.call(rbind, b.hat) - do.call(rbind, b)))
    }else if(diff.type == "abs.rel"){
      diffs <- abs(params.new-params)/(abs(params) + 1e-3)
      b.diff <- max(abs(do.call(rbind, b.hat) - do.call(rbind, b))/(abs(do.call(rbind, b)) + 1e-3))
    }
    diff <- max(diffs)
    # Message output (max relative diff)
    message("\nIteration ", iter + 1, " maximum difference: ", round(diff, 5))
    message("Largest change: ", names(params)[which(diffs==diff)])
    message("--> old: ", params[which(diffs==diff)], " new: ", params.new[which(diffs==diff)])
    message("Largest change in random effects: ", round(b.diff, 3))
    
    # Update ----
    params <- params.new
    D <- D.new; var.e <- var.e.new
    V <- lapply(mi, function(iii) {
      diag(x = rep(var.e, iii), ncol = sum(iii))
    })
    gamma <- gamma.new; gr <- rep(gamma, each = 2)
    eta <- eta.new
    beta <- beta.new; 
    b <- b.hat
    l0 <- l0.new; l0u <- l0u.new; l0i <- l0i.new; l0i.list <- as.list(l0i)
    iter <- iter + 1
    if(collect.hist) iter.hist = rbind(iter.hist, c(iter = iter, t(params)))
  }
  # Set up and return list ----
  coeffs <- list(beta = beta, var.e = var.e, D = D, gamma = gamma, eta = eta, hazard = cbind(ft, l0))
  
  rtn <- list(REs = do.call(rbind, b), coeffs = coeffs, 
              # Elapsed times //
              EMtime = round(sum(EM.time), 2),
              mvlme.time = mvlme.fit$elapsed.time,
              comp.time = round(proc.time()[3] - start.time, 2))
  
  if(post.process){
    message("\nStarting post-fit calculations...")
    pp.start <- proc.time()[3]
    
    # b at final parameter values.
    b <- mapply(function(b, Y, X, Z, V, mi, K, Delta, l0i, Fi, l0u, Fu, rvFi){
      ucminf::ucminf(b, ll, gradll,
                     Y, X, Z, V, D, sum(mi), K, Delta, l0i, Fi, l0u, Fu,
                     gamma, beta, eta, gr, rvFi, nK, q,
                     control = list(xtol = 1e-3, grtol = 1e-6))$par
    },b = b, Y = Y, X = X, Z = Z, V = V, mi = mi, K = K, Delta = Deltai.list, l0i = l0i.list,
    Fi = Fi.list, l0u = l0u, Fu = Fu, rvFi = rvFi.list, SIMPLIFY = F)
    bmat <- lapply(b, matrix, nc = 2, byr = T)
    
    # Covariance matrix for each subject at the MLEs for Omega and posterior mode bi.
    Sigmai <- mapply(function(b, Z, V, K, l0u, Fu){
      solve(-1 * sdll(b, Z, D, V, K, l0u, Fu, eta, gr, nK))
    }, b = b.hat, Z = Z, V = V, K = K, l0u = l0u, Fu = Fu, SIMPLIFY = F)
    S <- lapply(Sigmai, function(y) lapply(b.inds, function(x) y[x,x]))
    
    
    SEs <- hessian(coeffs, dmats, V, b, bmat, Sigmai, S, l0u, gh.nodes, n, q, nK)
    names(SEs) <- names(params)
    pp.end <- proc.time()[3]
    message("\nDone")
  }
  
  if(collect.hist) rtn$history <- iter.hist
  if(post.process){
    rtn$SEs <- SEs
    rtn$REs <- do.call(rbind, b)
    rtn$postprocess.time <- round(pp.end - pp.start, 2)
    rtn$comp.time = round(proc.time()[3] - start.time, 2)
  }
  rtn
}
