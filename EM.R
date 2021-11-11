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
  Di <- sv$Di
  l0 <- sv$l0; l0i <- sv$l0i; l0u <- sv$l0u
  Fi <- sv$Fi; Fu <- sv$Fu
  # List-versions of several data objects for use of mapply()
  Fi.list <- lapply(1:nrow(Fi), function(i) Fi[i, ])
  rvFi.list <- lapply(1:nrow(Fi), function(i) do.call(c, replicate(nK, Fi[i,], simplify = F)))
  Krep <- sapply(1:n, function(x){  # For updates to \eta
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
  l0.init <- inits.surv$l0.init
  
  # Extract all survival-related objects
  sv <- surv.mod(ph, data, l0.init)
  ft <- sv$ft; nev <- sv$nev
  surv.ids <- sv$surv.ids; surv.times <- sv$surv.times
  Di <- sv$Di
  l0 <- sv$l0; l0i <- sv$l0i; l0u <- sv$l0u
  Fi <- sv$Fi; Fu <- sv$Fu
  K <- getKi(data, ph)
  
  # Cast to parameter vector
  params <- c(vech(D), beta, var.e, gamma, eta)
  names(params) <- c(rep("D", length(vech(D))), 
                     names(beta), names(inits.long$var.e.init), names(gamma), paste0('eta_', names(eta)))
  
  # Collect history
  dmats <- list(Y=Y, X=X, Z=Z, mi=mi, Yk=Yk, Xk=Xk, Zk=Zk, K = K, sv=sv)
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
    },b = b, Y = Y, X = X, Z = Z, V = V, mi = mi, K = K, Delta = as.list(Di), l0i = as.list(l0i),
      Fi = Fi.list, l0u = l0u, Fu = Fu, rvFi = rvFi.list, SIMPLIFY = F)
    b.hat.split <- lapply(b.hat, function(y) lapply(b.inds, function(x) y[x]))
    
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
    #' Define \mu for survival submodel
    mu.surv <- mapply(function(K, Fu, b){
      rhs <- 0
      for(k in 1:nK) rhs <- rhs + gamma[k] * b[b.inds[[k]]]
      exp(K %*% eta + Fu %*% rhs)
    }, K = Krep, Fu = Fu, b = b.hat, SIMPLIFY = F)
    
    #' Define \tau for survival submodel
    tau <- mapply(function(Fu, S){
      out <- numeric(nrow(Fu))
      for(k in 1:nK) out <- out + diag(gamma[k]^2 * tcrossprod(Fu %*% S[[k]], Fu))
      out
    }, Fu = Fu, S = S, SIMPLIFY = F)
    
    tau.tilde <- mapply(function(Fu, S){
      mat <- matrix(0, nr = nK, nc = nrow(Fu))
      for(k in 1:nK) mat[k, ] <- diag(tcrossprod(Fu %*% S[[k]], Fu))
      mat
    }, Fu = Fu, S = S, SIMPLIFY = F)
    
    tau.surv <- lapply(tau, sqrt)
    tau2.surv <- lapply(tau, function(x){
      x <- x^(-0.5)
      if(any(is.infinite(x))) x[which(is.infinite(x))] <- 0   # Avoid NaN
      x
    })
    
    #' Set out Newton-Raphson items for update to (\gamma, \eta)
    #' S(\gamma) ----
    Sgamma <- mapply(function(Delta, tau.surv, mu.surv, l0u, Fu, Fi, b){
      t(Sgammacalc(gamma, Delta, tau.surv, mu.surv, l0u, Fu, Fi, w, v, b, nK, gh.nodes))
    }, Delta = as.list(Di), tau.surv = tau.surv, mu.surv = mu.surv, l0u = l0u,
    Fu = Fu, Fi = Fi.list, b = b.hat.split, SIMPLIFY = F)
    
    #' I(\gamma)
    Igamma <- mapply(function(tau.tilde, tau.surv, tau2.surv, mu.surv, Fu, l0u, b){
      gamma2Calc(gamma, tau.tilde, tau.surv, tau2.surv, mu.surv, w, v, Fu, l0u, b, nK, gh.nodes)
    }, tau.tilde = tau.tilde, tau.surv = tau.surv, tau2.surv = tau2.surv,
    mu.surv = mu.surv, Fu = Fu, l0u = l0u, b = b.hat.split, SIMPLIFY = F)
    
    #' S(\eta) ----
    Seta <- mapply(function(K, KK, Delta, l0u, mu.surv, tau.surv){
      cen <- Delta %*% K
      rhs <- c(0, 0)
      for(l in 1:gh.nodes) rhs <- rhs + w[l] * t(KK) %*% (l0u * (mu.surv * exp(tau.surv * v[l])))
      cen-t(rhs)
    }, K = K, KK = Krep, Delta = as.list(Di), l0u = l0u, mu.surv = mu.surv, tau.surv = tau.surv, SIMPLIFY = F)
    
    #' I(\eta) ----
    Ieta <- mapply(function(K, KK, tau.surv, mu.surv, l0u){
      Ietacalc(2, K, KK, tau.surv, mu.surv, l0u, w, v, gh.nodes)
    }, K = K, KK = Krep, tau.surv = tau.surv, mu.surv = mu.surv, l0u = l0u, SIMPLIFY = F)
    
    #' Second derivates of \gamma and \eta ('cross-terms') -----
    Igammaeta <- list() # for some reason I can't get mapply() to work with this 
    # Old
    for(i in 1:n){
      Igammaeta[[i]] <- Igammaetacalc(2, Krep[[i]], tau.surv[[i]], #tau.tilde[[i]],
                                      mu.surv[[i]], l0u[[i]], Fu[[i]], b.hat.split[[i]], gamma, w, v, nK, gh.nodes)
    }
    
    #' ##################
    #' M-step ========= #
    #' ##################
    
    #' D -----
    D.new <- Reduce('+', D.newi)/n
    
    #' \beta ----
    beta.new <- solve(Reduce('+', XtX)) %*% Reduce('+', beta.rhs) # NB this slightly faster than Reduce('+',.) on rhs
    #' var.e ----
    var.e.new <- colSums(do.call(rbind, Ee))/colSums(do.call(rbind, mi))      # NB this same speed as storing Ee directly as an array
    #' The baseline hazard, \lambda ----
    lambda <- lambdaUpdate(surv.times, ft, gamma, eta, K, S,
                           b.hat.split, n, w, v, gh.nodes, nK)                # (NB this not monitored for convergence)
    l0.new <- nev/rowSums(lambda)
    l0u.new <- lapply(l0u, function(x){
      ll <- length(x); l0.new[1:ll]
    })
    l0i.new <- c()
    l0i.new[which(Di == 0)] <- 0 
    l0i.new[which(Di == 1)] <- l0.new[match(Fi[which(Di==1), 2], ft)]
    
    #' (\gamma, \eta) ----
    # Set up score vector and information matrix
    Sge <- c(colSums(do.call(rbind, Sgamma)), colSums(do.call(rbind, Seta)))
    Imat <- as.matrix(Matrix::bdiag(Reduce('+', Igamma),
                                    Reduce('+', Ieta)))
    # Fill in off-block diagonal with Igammaeta
    eta.inds <- (nK+1):ncol(Imat)
    for(k in 1:nK){
      Imat[k, eta.inds] <- Imat[eta.inds, k] <- rowSums(do.call(cbind, lapply(Igammaeta, '[[', k)))
    }
    
    gamma.eta.new <- c(gamma, eta) + solve(Imat, Sge)
    
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
    l0 <- l0.new; l0u <- l0u.new; l0i <- l0i.new
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
    },b = b, Y = Y, X = X, Z = Z, V = V, mi = mi, K = K, Delta = as.list(Di), l0i = as.list(l0i),
    Fi = Fi.list, l0u = l0u, Fu = Fu, rvFi = rvFi.list, SIMPLIFY = F)
    
    # Covariance matrix for each subject at the MLEs for Omega and posterior mode bi.
    Sigmai <- mapply(function(b, Z, V, K, l0u, Fu){
      solve(-1 * sdll(b, Z, D, V, K, l0u, Fu, eta, gr, nK))
    }, b = b.hat, Z = Z, V = V, K = K, l0u = l0u, Fu = Fu, SIMPLIFY = F)
    
    SEs <- hessian(coeffs, dmats, V, b, Sigmai, l0i, l0u, gh.nodes, n, q, nK)
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
