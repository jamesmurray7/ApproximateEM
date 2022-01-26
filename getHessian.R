#' ####
#' Obtaining empirical information matrix for MLEs from an EM fit.
#' ----
#' Returns SEs for all parameters except for the baseline hazard.
#' ####

hessian <- function(Omega, data.mat, V, b, bmat, Sigmai, S, l0u, gh.nodes, n, q, nK){
  # Extract fitted values (MLEs) ----
  beta <- Omega$beta
  var.e <- Omega$var.e
  D <- Omega$D
  gamma <- Omega$gamma
  eta <- Omega$eta
  l0 <- Omega$hazard[, 2]

  # Extract data objects ----
  # Longitudinal //
  Z <- data.mat$Z; Zk <- data.mat$Zk
  X <- data.mat$X; Xk <- data.mat$Xk
  Y <- data.mat$Y; Yk <- data.mat$Yk
  mi <- data.mat$mi
  
  # Survival //
  K <- data.mat$K
  KK <- data.mat$KK
  Fi <- data.mat$Fi
  Fu <- data.mat$Fu
  Deltai <- data.mat$Deltai
  
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
      out <- 0.5 * tcrossprod(b) %*% (Dinv %*% delta.D[[i]] %*% Dinv)   # (Sigmai + tcrossprod(b))?
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
  svar.e <- list()
  for(i in 1:n){
    temp <- numeric(nK)
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
  Sge <- mapply(function(bmat, S, K, KK, Fu, Fi, l0u, Delta){
    Sgammaeta(c(gamma, eta), bmat, S, K, KK, Fu, Fi, l0u, Delta, w, v, 1e-4)
  }, bmat = bmat, S = S, K = K, KK = KK, Fu = Fu, Fi = Fi, 
  l0u = l0u, Delta = Deltai, SIMPLIFY = F)
  
  # Forming information -----------------------------------------------------
  si <- mapply(function(sD, Sb, Sv, Sge){
    c(sD, t(Sb), Sv, Sge)
  }, sD = sD, Sb = sbeta, Sv = svar.e, Sge = Sge)
  
  S <- rowSums(si)
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(si[, i]))) - S/n  # NB This RHS should = 0, however due to approximate nature of the approach, we leave this term in 
                                                                        # and note a small attenuation towards the null in SE calculation in the case it is absent.
  return(sqrt(diag(solve(I))))
}
