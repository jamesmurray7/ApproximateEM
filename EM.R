#' ####
#' Expectation Maximisation 
#' Assumes working directory is set to the parent directory
#'                    (i.e. ~/.../ApproximateEM)
#' ----
#' This is for K = 3, corresponding programs exist for K=5, and K=10, 
#'           (see the Kvariate folder).
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
#' 'cont' and 'bin' variables only (castData prepares data and ph in a suitable fashion)
#' 
#' The user can control the number of Gauss-Hermite quadrature nodes used (gh.nodes); the
#' maximum number of iterations to try before exiting (max.iter); the tolerance used
#' by the algorithm (tol); the convergence criterion to use (diff.type, must be 'abs' or
#' 'abs.rel'); if parameter estimates should be stored after each iteration (collect.hist);
#' two arguments control whether the optimisation step should no longer occur depending on 
#' iteration-to-iteration difference (b.tol) or a maximum number of optimisations done
#' (numOptim); nK denotes the number of covariates, but does NOT automate for K=5,10 shown
#' in the paper.
#' ####

source("./prerequisites.R")

em <- function(data, ph,
               gh.nodes = 3, collect.hist = T, max.iter = 200, 
               tol=0.01, diff.type = "abs.rel",
               b.tol = .1, numOptims = 4, nK = 3){
  # Set-up ------------------------------------------------------------------
  start.time <- proc.time()[3]
  diff <- 100; b.diff <- 100; iter <- 0
  uids <- unique(data$id); n <- length(uids)
  # Extract longitudinal-related objects
  X <- getXi(data, nK); Z <- getZi(data, nK); Y <- getYi(data, nK)
  mi <- getmi(data, nK)
  XtX <- lapply(X, crossprod)
  Zk <- splitZks(data, nK)
  Xk <- splitXks(data, nK)
  Yk <- splitYks(Y, mi, nK)
  inits.long <- suppressWarnings(Longit.inits(nK, data))
  inits.surv <- TimeVarCox(data, Ranefs(inits.long))
  
  # MVLME
  mvlme.fit <- mvlme(data, Y, X, Z, Yk, Xk, Zk, mi, inits.long, nK)
  
  # Extract initial conditions
  D <- mvlme.fit$D
  beta <- c(mvlme.fit$beta); beta.mat <- matrix(beta, nr = nK, byrow = T)
  var.e <- mvlme.fit$var.e
  V <- mvlme.fit$V
  mvlme.time <- mvlme.fit$elapsed.time
  gamma <- inits.surv$inits[3:length(inits.surv$inits)]
  eta <- inits.surv$inits[1:2]
  l0.init <- inits.surv$l0.init
  
  b0 <- matrix(NA, nr = length(uids), nc = 6)
  for(i in uids){
    b0[i, ] <- t(Z[[i]] %*% D) %*% solve(Z[[i]] %*% D %*% t(Z[[i]]) + V[[i]]) %*% (Y[[i]] - X[[i]] %*% beta)
  }
  
  sv <- surv.mod(ph, data, l0.init)
  
  # Extract all survival-related objects
  ft <- sv$ft; nev <- sv$nev
  surv.ids <- sv$surv.ids; surv.times <- sv$surv.times
  Di <- sv$Di
  l0 <- sv$l0; l0i <- sv$l0i; l0u <- sv$l0u
  Fi <- sv$Fi; Fu <- sv$Fu
  K <- getKi(data, ph)
  
  # Cast to parameter vector
  params <- c(vech(D), var.e, gamma, beta, eta)
  names(params) <- c(rep("D", length(vech(D))), "var.e_1", "var.e_2", "var.e_3",
                     paste0("gamma_", 1:3),
                     paste0("beta1_", 0:3), paste0("beta2_", 0:3), paste0("beta3_", 0:3),
                     paste0("surv.", names(eta)))
  
  # Collect history
  if(collect.hist) iter.hist = data.frame(iter = iter, t(params))
  
  # Gaussian Quadrature -----------------------------------------------------
  gh <- statmod::gauss.quad.prob(gh.nodes, "normal")
  v <- gh$nodes
  w <- gh$weights
  message("Starting EM Algorithm")
  EM.time <- c()
  # EM ----------------------------------------------------------------------
  while(diff > tol & iter < max.iter){
    # E-step ----
    Eresid <- Sgamma <- matrix(NA, nrow = length(uids), nc = 3)
    Seta <- Igamma1eta <- Igamma2eta <- Igamma3eta <- matrix(NA, nrow = length(uids), nc = 2)
    D.newi <- Sigmai.store <- beta.part <- Ieta <- Igamma <-  list()
    bi.mat <- matrix(NA, nrow = length(uids), nc = 6)
    lambda.store <- matrix(0, nrow = length(ft), nc = length(uids))
    # Begin loop over subjects
    p1 <- proc.time()[3]
    for(i in uids){
      if(iter <= numOptims || b.diff > b.tol){
        bi <- ucminf::ucminf(b0[i,], bllC, gradC,
                             Y[[i]], X[[i]], Z[[i]], V[[i]], D, nrow(Z[[i]]),
                             K[[i]], Di[i], l0i[i], Fi[i, ], l0u[[i]],  Fu[[i]], gamma, beta, eta,
                             rep(gamma, each = 2), repVec(Fi[i, ]),
                             hessian = 0, control = list(xtol = 1e-3, grtol = 1e-6))$par
      }else{
        bi <- b0[i, ]
      }
      
      Sigmai <- solve(-1 * d2b.ll(K[[i]], Z[[i]], D, V[[i]], l0u[[i]], Fu[[i]], gamma, eta, bi))
      
      # Collect the random effects
      bi.mat[i,] <- bi; b1 <- bi[1:2]; b2 <- bi[3:4]; b3 <- bi[5:6]
      bb <- rbind(b1, b2, b3)
      # And the variance
      Sigmai.store[[i]] <- Sigmai
      S <- list()
      # Covariance matrices on the first, second and third REs
      S[[1]] <- Sigmai[1:2,1:2]; S[[2]] <- Sigmai[3:4, 3:4]; S[[3]] <- Sigmai[5:6, 5:6]
      
      # Step to update D
      D.newi[[i]] <- Sigmai + tcrossprod(bi)
      
      # Step to update \sigma_\epsilon
      Eresidi <- matrix(NA, nr = gh.nodes, nc = 3)
      for(L in 1:3){
        for(k in 1:gh.nodes){
          mu.long.K <- Xk[[i]][[L]] %*% beta.mat[L,] + Zk[[i]][[L]] %*% bb[L,]
          tau.long.K <- sqrt(diag(Zk[[i]][[L]] %*% S[[L]] %*% t(Zk[[i]][[L]])))
          Eresidi[k, L] <- w[k] * crossprod(Yk[[i]][,L] - mu.long.K - tau.long.K * v[k])
        }
      }
      Eresid[i,] <- colSums(Eresidi)
      
      # Step to update \beta
      beta.rhs <- 0
      mu.long <- Z[[i]] %*% bi
      tau.long <- sqrt(diag(Z[[i]] %*% Sigmai %*% t(Z[[i]])))
      for(k in 1:gh.nodes) beta.rhs <- beta.rhs + w[k] * tau.long * v[k]
      beta.part[[i]] <- crossprod(X[[i]], Y[[i]] - mu.long - beta.rhs)
      
      # Update for \lambda_0 ----
      # Populate lambda.store matrix
      st <- surv.times[[i]]
      if(length(st) > 0){
        for(u in st){
          temp <- 0
          Fst <- cbind(1, ft[u])
          mu.surv <- exp(K[[i]] %*% eta + Fst %*% (gamma[1] * b1 + gamma[2] * b2 + gamma[3] * b3))
          taust <- sqrt(diag(
            gamma[1]^2 * tcrossprod(Fst %*% S[[1]], Fst) +
            gamma[2]^2 * tcrossprod(Fst %*% S[[2]], Fst) +
            gamma[3]^2 * tcrossprod(Fst %*% S[[3]], Fst)
          ))
          for(k in 1:gh.nodes) temp <- temp + w[k] * mu.surv * exp(v[k] * taust)
          lambda.store[u,i] <- temp
        }
      }
      
      # Survival parameter updates ----
      tau.tilde <- matrix(0, nr = 3, nc = nrow(Fu[[i]]))
      tau.tilde[1,] <- diag(tcrossprod(Fu[[i]] %*% S[[1]], Fu[[i]]))
      tau.tilde[2,] <- diag(tcrossprod(Fu[[i]] %*% S[[2]], Fu[[i]]))
      tau.tilde[3,] <- diag(tcrossprod(Fu[[i]] %*% S[[3]], Fu[[i]]))
      mu.surv <- exp(K[[i]] %*% eta) %x% exp(Fu[[i]] %*% (gamma[1] * b1 + gamma[2] * b2 + gamma[3] * b3))
      tau.surv <- sqrt(diag(gamma[1]^2 * tcrossprod(Fu[[i]] %*% S[[1]], Fu[[i]])) + 
                         diag(gamma[2]^2 * tcrossprod(Fu[[i]] %*% S[[2]], Fu[[i]])) + 
                         diag(gamma[3]^2 * tcrossprod(Fu[[i]] %*% S[[3]], Fu[[i]])))
      tau2.surv <- (diag(gamma[1]^2 * tcrossprod(Fu[[i]] %*% S[[1]], Fu[[i]])) + 
                      diag(gamma[2]^2 * tcrossprod(Fu[[i]] %*% S[[2]], Fu[[i]])) + 
                      diag(gamma[3]^2 * tcrossprod(Fu[[i]] %*% S[[3]], Fu[[i]])))^(-0.5)
      
      if(any(is.infinite(tau2.surv))){ # avoid NaN
        tau2.surv[which(is.infinite(tau2.surv))] <- 0
      }
      # S(\gamma) ----
      Sgamma.store <- matrix(NA, nr = gh.nodes, nc = 3)
      Sgammai <- c()
      for(L in 1:3){
        for(k in 1:gh.nodes){
          xi <- l0u[[i]] * (mu.surv * exp(v[k] * tau.surv))
          Sgamma.store[k, L] <- w[k] * crossprod(xi, Fu[[i]] %*% bb[L,]) + 
            gamma[L] * v[k] * w[k] * crossprod(xi * tau.surv, xi)
        }
        cen <- Di[i] * Fi[i,] %*% bb[L,]
        Sgammai[L] <- cen - sum(Sgamma.store[, L])
      }
      Sgamma[i,] <- Sgammai
      
      Igamma[[i]] <- gamma2Calc(gamma, tau.tilde, tau.surv, tau2.surv, mu.surv, w, v, Fu[[i]], l0u[[i]], bb, L = 3, gh = gh.nodes)
      
      # S(\eta) ----
      KK <- apply(K[[i]], 2, rep, nrow(Fu[[i]]))
      if("numeric" %in% class(KK)) KK <- t(as.matrix(KK))
      cen <- Di[i] %*% K[[i]]
      Seta.store <- matrix(NA, nrow = gh.nodes, nc = 2)
      for(k in 1:gh.nodes) Seta.store[k,] <- w[k] * t(KK) %*% (l0u[[i]] * (mu.surv * exp(tau.surv * v[k])))
      Seta[i,] <- cen - colSums(Seta.store)
      
      # I(\eta) ----
      Ieta.store <- list()
      for(k in 1:gh.nodes){
        if(length(mu.surv) == 1){
          diagvec <- diag(l0u[[i]] * (mu.surv * exp(tau.surv * v[k])))
        }else{
          diagvec <- diag(an(l0u[[i]] * (mu.surv * exp(tau.surv * v[k]))))
        }
        Ieta.store[[k]] <- w[k] * crossprod(diagvec %*% KK, KK)
      }
      
      Igamma1eta.store <- Igamma2eta.store <- Igamma3eta.store <- matrix(NA, nr = 2, nc = gh.nodes)
      for(k in 1:gh.nodes){
        xi <- l0u[[i]] * (mu.surv * exp(v[k] * tau.surv))
        Igamma1eta.store[,k] <- w[k] * crossprod(KK, (Fu[[i]] %*% bb[1,]) * xi) + 
          2 * gamma[1] * w[k] * v[k] * crossprod(KK, xi * tau.surv * xi)
        Igamma2eta.store[,k] <- w[k] * crossprod(KK, (Fu[[i]] %*% bb[2,]) * xi) + 
          2 * gamma[2] * w[k] * v[k] * crossprod(KK, xi * tau.surv * xi)
        Igamma3eta.store[,k] <- w[k] * crossprod(KK, (Fu[[i]] %*% bb[3,]) * xi) + 
          2 * gamma[3] * w[k] * v[k] * crossprod(KK, xi * tau.surv * xi)
      }
      
      Ieta[[i]] <- Reduce('+', Ieta.store)
      Igamma1eta[i,] <- rowSums(Igamma1eta.store)
      Igamma2eta[i,] <- rowSums(Igamma2eta.store)
      Igamma3eta[i,] <- rowSums(Igamma3eta.store)
    }
    
    # M-step ----
    
    # D
    D.new <- Reduce('+', D.newi)/n
    # \beta and \epsilon_k
    beta.new <- solve(Reduce('+', XtX)) %*% Reduce('+', beta.part)
    var.e.new <- colSums(Eresid)/colSums(do.call(rbind, mi))
    # \lambda_0
    l0.new <- nev/rowSums(lambda.store)
    l0u.new <- lapply(l0u, function(x){
      ll <- length(x); l0.new[1:ll]
    })
    l0i.new <- c()
    l0i.new[which(Di == 0)] <- 0 
    l0i.new[which(Di == 1)] <- l0.new[match(Fi[which(Di==1), 2], ft)]
    
    # \gamma and \eta ----
    # Setting up score and information
    Sge <- c(colSums(Sgamma), colSums(Seta))
    Imat <- as.matrix(Matrix::bdiag(Reduce('+', Igamma), Reduce('+', Ieta)))
    Imat[1,4:5] <- Imat[4:5,1] <- colSums(Igamma1eta)
    Imat[2,4:5] <- Imat[4:5,2] <- colSums(Igamma2eta)
    Imat[3,4:5] <- Imat[4:5,3] <- colSums(Igamma3eta)
    
    gamma.eta.new <- c(gamma, an(eta)) + solve(Imat, Sge)
    gamma.new <- gamma.eta.new[1:3]
    eta.new <- gamma.eta.new[4:5]
    
    EM.time[iter + 1] <- proc.time()[3] - p1 # M step finishes here ---
    
    # Collect, take differences and print ----
    params.new <- c(vech(D.new), var.e.new, gamma.new, beta.new, eta.new)
    # Take differences (user input)
    if(diff.type == "abs"){
      diffs <- abs(params.new-params)
    }else if(diff.type == "abs.rel"){
      diffs <- abs(params.new-params)/abs(params)
    }
    diff <- max(diffs)
    # Take differences (user input)
    if(diff.type == "abs"){
      diffs <- abs(params.new-params)
      b.diff <- max(abs(bi.mat - b0))
    }else if(diff.type == "abs.rel"){
      diffs <- abs(params.new-params)/(abs(params) + 1e-3)
      b.diff <- max(abs(bi.mat - b0)/(abs(b0) + 1e-3))
    }
    diff <- max(diffs)
    # Print
    # print(sapply(params, round, 4))
    message("\nIteration ", iter + 1, " maximum difference: ", round(diff, 5))
    message("Largest change: ", names(params)[which(diffs==diff)])
    message("--> old: ", params[which(diffs==diff)], " new: ", params.new[which(diffs==diff)])
    
    message("Largest change in random effects: ", round(b.diff, 3))
    
    # Update ----
    names(params.new) <- names(params); params <- params.new
    D <- D.new; var.e <- var.e.new
    V <- lapply(mi, function(iii) {
      diag(x = rep(var.e, iii), ncol = sum(iii))
    })
    gamma <- gamma.new
    eta <- eta.new
    beta <- beta.new; beta.mat <- rbind(beta[1:4], beta[5:8], beta[9:12])
    b0 <- bi.mat
    l0 <- l0.new; l0u <- l0u.new; l0i <- l0i.new
    iter <- iter + 1
    if(collect.hist) iter.hist = rbind(iter.hist, c(iter = iter, t(params)))
  }
  # Set up and return list ----
  rtn <- list(REs = bi.mat, var.e = var.e, D = D, gamma = gamma, hazard = cbind(ft, l0),
              beta = beta, eta = eta, EMtime = round(sum(EM.time), 2),
              mvlme.time = mvlme.fit$elapsed.time,
              comp.time = round(proc.time()[3] - start.time, 2))
  if(collect.hist) rtn$history <- iter.hist
  rtn
}

