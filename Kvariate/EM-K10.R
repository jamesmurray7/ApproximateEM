#' ####
#' Expectation Maximisation (K=10).
#' Assumes working directory is set to the parent directory
#'                    (i.e. ~/.../ApproximateEM)
#' ----
#' This is for K = 10, 
#' All arguments and calls are as per the EM.R shown for K=3 (in parent directory),
#' see program header there for details.
#' This begins (admittedly inefficiently) by re-defining functions from prerequisites.R.
#' ####

source("./prerequisites.R")
# Overwrite the loglik for K=3.
sourceCpp("./source/ll10.cpp")

# Overwrite second derivative - needed for \Sigma (This should just be an updated function to
# include K!).
d2b.ll <- function(K, Z, D, V, l0u, Fu, g, eta, bi){
  gr <- rep(g, each = 2)
  b <- bi
  diag.g <- diag(gr)
  if(nrow(Fu) == 1){
    surv.part <- -diag.g %*% repCols(Fu, 10) %*% (l0u * exp(K %*% eta + repCols(Fu, 10) %*% (gr * b))) %*% repCols(Fu, 10) %*% diag.g
  }else{
    diag.kern <- diag(l0u * an(exp(K %*% eta) %x% exp(repCols(Fu, 10) %*% (gr * b))))
    surv.part <- -diag.g %*% crossprod(repCols(Fu, 10), diag.kern) %*% repCols(Fu, 10) %*% diag.g
  }
  -crossprod(Z, solve(V) %*% Z)- solve(D) + surv.part
}

# The EM Algorithm --------------------------------------------------------

em <- function(data, ph,
               gh.nodes = 3, collect.hist = T, max.iter = 200, 
               tol=0.01, diff.type = "abs.rel",
               b.tol = .1, numOptims = 4, nK = 10){
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
  
  b0 <- matrix(NA, nr = length(uids), nc = nK * 2)
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
  names(params) <- c(rep("D", length(vech(D))), 
                     paste0("var.e_", 1:nK),
                     paste0("gamma_", 1:nK),
                     paste0(rep(paste0("beta", 1:nK), each = 4), 0:3),
                     paste0("surv.", names(eta)))
  params
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
    Eresid <- Sgamma <- matrix(NA, nrow = length(uids), nc = 10)
    Seta <- Igamma1eta <- Igamma2eta <- Igamma3eta <- Igamma4eta <- Igamma5eta <- matrix(NA, nrow = length(uids), nc = 2)
    Igamma6eta <- Igamma7eta <- Igamma8eta <- Igamma9eta <- Igamma10eta <- Igamma5eta
    D.newi <- Sigmai.store <- beta.part <- Ieta <- Igamma <- list()
    bi.mat <- matrix(NA, nrow = length(uids), nc = 20)
    lambda.store <- matrix(0, nrow = length(ft), nc = length(uids))
    # Begin loop over subjects
    p1 <- proc.time()[3]
    for(i in uids){
      if(iter <= numOptims || b.diff > b.tol){
        bi <- ucminf::ucminf(b0[i,], bllC, gradC,
                             Y[[i]], X[[i]], Z[[i]], V[[i]], D, nrow(Z[[i]]),
                             K[[i]], Di[i], l0i[i], Fi[i, ], l0u[[i]],  Fu[[i]], gamma, beta, eta,
                             rep(gamma, each = 2), repVec(Fi[i, ], 10),
                             hessian = 0, control = list(xtol = 1e-3, grtol = 1e-6))$par
      }else{
        bi <- b0[i, ]
      }
      
      Sigmai <- solve(-1 * d2b.ll(K[[i]], Z[[i]], D, V[[i]], l0u[[i]], Fu[[i]], gamma, eta, bi))
      
      # Collect the random effects
      bi.mat[i,] <- bi
      b1 <- bi[1:2]; b2 <- bi[3:4]; b3 <- bi[5:6]; b4 <- bi[7:8]; b5 <- bi[9:10]
      b6 <- bi[11:12]; b7 <- bi[13:14]; b8 <- bi[15:16]; b9 <- bi[17:18]; b10 <- bi[19:20]
      bb <- rbind(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
      # And the variance
      Sigmai.store[[i]] <- Sigmai
      # Covariance matrices on the first, second and third REs
      S <- list()
      # This automates, but is about 3x slower (microsecs).
      # S <- lapply(split(seq(nK * 2), rep(1:nK, each = 2)), function(x) Sigmai[x,x])
      S[[1]] <- Sigmai[1:2,1:2]; S[[2]] <- Sigmai[3:4, 3:4]; S[[3]] <- Sigmai[5:6, 5:6]
      S[[4]] <- Sigmai[7:8,7:8]; S[[5]] <- Sigmai[9:10,9:10]; S[[6]] <- Sigmai[11:12, 11:12]
      S[[7]] <- Sigmai[13:14, 13:14]; S[[8]] <- Sigmai[15:16,15:16]; S[[9]] <- Sigmai[17:18, 17:18]
      S[[10]] <- Sigmai[19:20, 19:20]
      
      # Step to update D
      D.newi[[i]] <- Sigmai + tcrossprod(bi)
      
      # Step to update \sigma_\epsilon
      Eresidi <- matrix(NA, nr = gh.nodes, nc = 10)
      for(L in 1:10){
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
          mu.surv <- exp(K[[i]] %*% eta + Fst %*% (gamma[1] * b1 + gamma[2] * b2 + gamma[3] * b3 + 
                                                     gamma[4] * b4 + gamma[5] * b5 + gamma[6] * b6 + 
                                                     gamma[7] * b7 + gamma[8] * b8 + gamma[9] * b9 + 
                                                     gamma[10] * b10))
          taust <- sqrt(diag(
            gamma[1]^2 * tcrossprod(Fst %*% S[[1]], Fst) + 
              gamma[2]^2 * tcrossprod(Fst %*% S[[2]], Fst) + 
              gamma[3]^2 * tcrossprod(Fst %*% S[[3]], Fst) + 
              gamma[4]^2 * tcrossprod(Fst %*% S[[4]], Fst) + 
              gamma[5]^2 * tcrossprod(Fst %*% S[[5]], Fst) + 
              gamma[6]^2 * tcrossprod(Fst %*% S[[6]], Fst) + 
              gamma[7]^2 * tcrossprod(Fst %*% S[[7]], Fst) + 
              gamma[8]^2 * tcrossprod(Fst %*% S[[8]], Fst) + 
              gamma[9]^2 * tcrossprod(Fst %*% S[[9]], Fst) + 
              gamma[10]^2 * tcrossprod(Fst %*% S[[10]], Fst)
          ))
          for(k in 1:gh.nodes) temp <- temp + w[k] * mu.surv * exp(v[k] * taust)
          lambda.store[u,i] <- temp
        }
      }
      
      # Survival parameter updates ----
      tau.tilde <- matrix(0, nr = 10, nc = nrow(Fu[[i]]))
      for(L in 1:10){
        tau.tilde[L, ] <- diag(tcrossprod(Fu[[i]] %*% S[[L]], Fu[[i]]))
      }
      mu.surv <- exp(K[[i]] %*% eta) %x% exp(Fu[[i]] %*% (gamma[1] * b1 + gamma[2] * b2 + gamma[3] * b3 + 
                                                            gamma[4] * b4 + gamma[5] * b5 + gamma[6] * b6 + 
                                                            gamma[7] * b7 + gamma[8] * b8 + gamma[9] * b9 + 
                                                            gamma[10] * b10))
      tau <- diag(gamma[1]^2 * tcrossprod(Fu[[i]] %*% S[[1]], Fu[[i]])) + 
        diag(gamma[2]^2 * tcrossprod(Fu[[i]] %*% S[[2]], Fu[[i]])) + 
        diag(gamma[3]^2 * tcrossprod(Fu[[i]] %*% S[[3]], Fu[[i]])) + 
        diag(gamma[4]^2 * tcrossprod(Fu[[i]] %*% S[[4]], Fu[[i]])) +
        diag(gamma[5]^2 * tcrossprod(Fu[[i]] %*% S[[5]], Fu[[i]])) + 
        diag(gamma[6]^2 * tcrossprod(Fu[[i]] %*% S[[6]], Fu[[i]])) + 
        diag(gamma[7]^2 * tcrossprod(Fu[[i]] %*% S[[7]], Fu[[i]])) + 
        diag(gamma[8]^2 * tcrossprod(Fu[[i]] %*% S[[8]], Fu[[i]])) + 
        diag(gamma[9]^2 * tcrossprod(Fu[[i]] %*% S[[9]], Fu[[i]])) +
        diag(gamma[10]^2 * tcrossprod(Fu[[i]] %*% S[[10]], Fu[[i]]))  
      
      tau.surv <- sqrt(tau)
      tau2.surv <- tau^(-0.5)
      
      if(any(is.infinite(tau2.surv))){ # avoid NaN
        tau2.surv[which(is.infinite(tau2.surv))] <- 0
      }
      
      # S(\gamma) ----
      Sgamma.store <- matrix(NA, nr = gh.nodes, nc = 10)
      Sgammai <- c()
      for(L in 1:10){
        for(k in 1:gh.nodes){
          xi <- l0u[[i]] * (mu.surv * exp(v[k] * tau.surv))
          Sgamma.store[k, L] <- w[k] * crossprod(xi, Fu[[i]] %*% bb[L,]) + 
            gamma[L] * v[k] * w[k] * crossprod(xi * tau.surv, xi)
        }
        cen <- Di[i] * Fi[i,] %*% bb[L,]
        Sgammai[L] <- cen - sum(Sgamma.store[, L])
      }
      Sgamma[i,] <- Sgammai
      
      Igamma[[i]] <- gamma2Calc(gamma, tau.tilde, tau.surv, tau2.surv, mu.surv, w, v, Fu[[i]], l0u[[i]], bb, 10, 3)
      
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
      
      Igammaeta.store <- list()
      for(L in 1:10){
        Igammaeta.store[[L]] <- c(0, 0)
        for(k in 1:3){
          xi <- l0u[[i]] * (mu.surv * exp(v[k] * tau.surv))
          Igammaeta.store[[L]] <- Igammaeta.store[[L]] + w[k] * crossprod(KK, (Fu[[i]] %*% bb[L,]) * xi) + 
            2 * gamma[L] * w[k] * v[k] * crossprod(KK, xi * tau.surv * xi)
        }
      }
      
      Ieta[[i]] <- Reduce('+', Ieta.store)
      Igamma1eta[i,] <- Igammaeta.store[[1]]; Igamma2eta[i,] <- Igammaeta.store[[2]]
      Igamma3eta[i,] <- Igammaeta.store[[3]]; Igamma4eta[i,] <- Igammaeta.store[[4]]
      Igamma5eta[i,] <- Igammaeta.store[[5]]; Igamma6eta[i,] <- Igammaeta.store[[6]]
      Igamma7eta[i,] <- Igammaeta.store[[7]]; Igamma8eta[i,] <- Igammaeta.store[[8]]
      Igamma9eta[i,] <- Igammaeta.store[[9]]; Igamma10eta[i,] <- Igammaeta.store[[10]]
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
    l0i.new[which(Di == 1)] <- l0.new[match(unique(data[data$status == 1,]$survtime), ft)]
    
    # \gamma and \eta ----
    # Setting up score and information
    Sge <- c(colSums(Sgamma), colSums(Seta))
    Imat <- as.matrix(Matrix::bdiag(Reduce('+', Igamma), Reduce('+', Ieta)))
    Imat[1,11:12] <- Imat[11:12,1] <- colSums(Igamma1eta)
    Imat[2,11:12] <- Imat[11:12,2] <- colSums(Igamma2eta)
    Imat[3,11:12] <- Imat[11:12,3] <- colSums(Igamma3eta)
    Imat[4,11:12] <- Imat[11:12,4] <- colSums(Igamma4eta)
    Imat[5,11:12] <- Imat[11:12,5] <- colSums(Igamma5eta)
    Imat[6,11:12] <- Imat[11:12,6] <- colSums(Igamma6eta)
    Imat[7,11:12] <- Imat[11:12,7] <- colSums(Igamma7eta)
    Imat[8,11:12] <- Imat[11:12,8] <- colSums(Igamma8eta)
    Imat[9,11:12] <- Imat[11:12,9] <- colSums(Igamma9eta)
    Imat[10,11:12] <- Imat[11:12,10] <- colSums(Igamma10eta)
    
    gamma.eta.new <- c(gamma, an(eta)) + solve(Imat, Sge)
    gamma.new <- gamma.eta.new[1:10]
    eta.new <- gamma.eta.new[11:12]
    
    EM.time[iter + 1] <- proc.time()[3] - p1 # M step finishes here ---
    
    # Collect, take differences and print ----
    params.new <- c(vech(D.new), var.e.new, gamma.new, beta.new, eta.new)
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
    beta <- beta.new; beta.mat <- matrix(beta, nr = nK, byrow = T)
    b0 <- bi.mat
    l0 <- l0.new; l0u <- l0u.new; l0i <- l0i.new
    iter <- iter + 1
    if(collect.hist) iter.hist = rbind(iter.hist, c(iter = iter, t(params)))
  }
  # Set up and return list ----
  rtn <- list(REs = bi.mat, var.e = var.e, D = D, gamma = gamma, hazard = cbind(ft, l0),
              beta = beta, eta = eta, elapsed.time = round(sum(EM.time), 2))
  if(collect.hist) rtn$history <- iter.hist
  rtn
}
