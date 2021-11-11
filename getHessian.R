#' ####
#' Obtaining empirical information matrix for MLEs from an EM fit.
#' ----
#' Returns SEs for all parameters except for the baseline hazard.
#' ####

hessian <- function(Omega, data.mat, V, b, S, l0i, l0u, gh.nodes, n, q, nK){
  # Extract fitted values (MLEs) ----
  beta <- Omega$beta
  var.e <- Omega$var.e
  D <- Omega$D
  gamma <- Omega$gamma
  eta <- Omega$eta
  l0 <- Omega$hazard[, 2]
  Sigmai <- S
  
  # Extract data objects ----
  # Longitudinal //
  Z <- data.mat$Z; Zk <- data.mat$Zk
  X <- data.mat$X; Xk <- data.mat$Xk
  Y <- data.mat$Y; Yk <- data.mat$Yk
  mi <- data.mat$mi
  # Survival //
  K <- data.mat$K
  sv <- data.mat$sv
  Fi <- sv$Fi; Fu <- sv$Fu; Deltai <- sv$Di; ft <- sv$ft
  
  gh <- statmod::gauss.quad.prob(gh.nodes, 'normal') 
  w <- gh$weights; v <- gh$nodes
  
  # Scores ----------------------------------------------------------------
  # S(D) ----
  # All credit to Graeme Hickey
  # https://github.com/graemeleehickey/joineRML/blob/master/R/hessian.R (lines 66 - 92).
  
  Dinv <- solve(D)
  vech.indices <- which(lower.tri(D, diag = T), arr.ind = T)
  dimnames(vech.indices) <- NULL
  delta.D <- lapply(1:nrow(vech.indices), function(d){
    out <- matrix(0, nrow(D), ncol(D))
    ind <- vech.indices[d, 2:1]
    out[ind[1], ind[2]] <- out[ind[2], ind[1]] <- 1 # dD/dvech(d)_i
    out
  })
  
  lhs <- sapply(delta.D, function(d) {
    -0.5 * sum(diag(Dinv %*% d))
  })
  
  sDi <- function(i) {
    mapply(function(b) {
      out <- 0.5 * tcrossprod(b) %*% (Dinv %*% delta.D[[i]] %*% Dinv)
      lhs[i] + sum(diag(out))
    },
    b = b,
    SIMPLIFY = T)
  }
  
  sD <- sapply(1:nrow(vech.indices), sDi)
  sD <- lapply(1:nrow(sD), function(x) sD[x, ]) # Cast to list
  
  # Longitudinal ----
  # Define indices for beta and b
  b.inds <- split(seq(q), cut(seq_along(seq(q)), nK, labels = F))
  beta.inds <- split(seq(length(beta)), cut(seq_along(seq(length(beta))), nK, labels = F))
  
  # S(\beta) -----
  tau.long <- mapply(function(Z, S){
    sqrt(diag(Z %*% S %*% t(Z)))
  },
  Z = Z, S = Sigmai, SIMPLIFY = F
  )
  
  sbeta <- mapply(function(X, V, Y, Z, b, tau){
    rhs <- 0
    for(i in 1:3) rhs <- rhs + w[i] * v[i] * tau
    crossprod(X, solve(V) %*% (Y - X %*% beta - Z %*% b - rhs))
  },
  X = X, V = V, Y = Y, Z = Z, b = b, tau = tau.long, SIMPLIFY = F)
  
  # S(var.e) ----
  svar.e <- list()#matrix(NA, nr = n, nc = nK)
  for(i in 1:n){
    temp <- numeric(3)
    for(k in 1:nK){
      tauK <- sqrt(diag(Zk[[i]][[k]] %*% Sigmai[[i]][b.inds[[k]], b.inds[[k]]] %*% t(Zk[[i]][[k]])))
      rhs <- 0
      for(l in 1:gh.nodes){
        rhs <- rhs + w[l] * crossprod(Yk[[i]][, k] - Xk[[i]][[k]] %*% beta[beta.inds[[k]]] - 
                                        Zk[[i]][[k]] %*% b[[i]][b.inds[[k]]] - v[l] * tauK)
        
      } 
      temp[k] <- -mi[[i]][k]/(2 * var.e[k]) + 1/(2 * var.e[k]^2) * rhs
    }
    svar.e[[i]] <- temp
  }
  
  # Survival parameters ----
  # Define mu, tau and sum(gamma_k*b_{ik}) for use in Scores for gamma and eta
  tau.surv <- mapply(function(b, S, Fu){
    out <- 0
    for(k in 1:nK){
      out <- out + diag(gamma[k]^2 * tcrossprod(Fu %*% S[b.inds[[k]], b.inds[[k]]], Fu))
    }
    sqrt(out)
  }, b = b, S = Sigmai, Fu = Fu, SIMPLIFY = F)
  
  sumgammabk <- mapply(function(b){
    out <- 0
    for(k in 1:nK){
      out <- out + gamma[k] * b[b.inds[[k]]]
    }
    out
  }, b = b, SIMPLIFY = F)
  
  mu.surv <- mapply(function(K, Fu, gb){
    exp(K %*% eta) %x% exp(Fu %*% gb)
  }, K = K, Fu = Fu, gb = sumgammabk, SIMPLIFY = F)
  
  # S(gamma) -----
  sgamma <- list()
  for(i in 1:n){
    out <- matrix(rep(0, nK), nr = 1, nc = nK)
    for(k in 1:nK){
      cenK <- Deltai[i] * Fi[i, ] %*% b[[i]][b.inds[[k]]]
      rhs <- 0
      for(l in 1:3){
        xi <- getxi(tau.surv[[i]], mu.surv[[i]], v[l], l0u[[i]])
        rhs <- rhs + w[l] * crossprod(xi, Fu[[i]] %*% b[[i]][b.inds[[k]]]) + 
                      gamma[k] * w[l] * v[l] * crossprod(xi * tau.surv[[i]], xi)
      }
      out[1, k] <- cenK - rhs
    }
    sgamma[[i]] <- out
  }

  # S(eta) ----
  KK <- mapply(function(K, Fu){
    KK <- apply(K, 2, rep, nrow(Fu))
    if("numeric" %in% class(KK)) KK <- t(as.matrix(KK))
    KK
  }, K = K, Fu = Fu, SIMPLIFY = F)
  
  seta <- list()#matrix(NA, nr = n, nc = 2)
  for(i in 1:n){
    rhs <- 0
    for(l in 1:3) rhs <- rhs + w[l] * l0u[[i]] * exp(KK[[i]] %*% eta + Fu[[i]] %*% sumgammabk[[i]] + tau.surv[[i]] * v[l])
    seta[[i]] <- Deltai[i] %*% K[[i]] - t(crossprod(KK[[i]], rhs))
  }
  
  # Forming information -----------------------------------------------------
  Si <- mapply(function(sD, Sb, Sv, Sg, Se){
    c(sD, t(Sb), Sv, Sg, Se)
  }, sD = sD, Sb = sbeta, Sv = svar.e, Sg = sgamma, Se = seta, SIMPLIFY = F)
  
  lhs <- Reduce('+', lapply(Si, tcrossprod))
  rhs <- tcrossprod(colSums(do.call(rbind, Si)))
  I <- lhs - rhs/n    # NB This RHS should = 0, however due to approximate nature of the approach, we leave this term in 
                      # and note a small attenuation towards the null in SE calculation in the case it is absent.
  H <- solve(I)
  
  return(sqrt(diag(H)))
}
