###############################################################################
# update IBP parameter alpha
###############################################################################

update_alpha = function(K, H_D, e, f){
  # K  : current number of active components
  # H_D: D-th harmonic number (this should be calculated only once!)
  # e  : shape parameter of prior for alpha
  # f  : rate parameter of prior for alpha
  
  rgamma(n = 1L, shape = e+K, rate = f+H_D)
}

###############################################################################
# update for factors' precisions lambda
###############################################################################

update_lambda = function(G, Z, c, d){
  # Z: (truncated) IBP matrix, of size D \times K
  # G: mixing matrix, of size D \times K_+
  # c: shape parameter of prior for each lambda_k
  # d: rate parameter of prior for each lambda_k
  
  K = ncol(Z)
  
  # get matrix of parameters size K \times 2 (shape and rate)
  l_pars = matrix(NA_real_, nrow = K, ncol = 2L)
  colnames(l_pars) = c("shape", "rate")
  
  # get vector of shape parameters (length K)
  l_pars[,"shape"] = c + 0.5*colSums(Z)
  
  # get vector of rate parameters (length K)
  l_pars[,"rate"] = d + 0.5*colSums(G*G)
  
  # draw lambda
  lambda = apply(l_pars, 1L, function(p){
    rgamma(n=1L, shape = p["shape"], rate = p["rate"])
  })
  
  return(lambda)
}

###############################################################################
# update for lambda prior parameter d
###############################################################################

update_d = function(lambda, c, c_0, d_0){
  # lambda: vector with precisions lambda_k for sampling g_dk, length K
  # c: shape parameter of prior for each lambda_k
  # c_0: shape parameter of prior for d
  # d_0: rate parameter of prior for d
  
  rgamma(n=1L, shape = c_0 + c*length(lambda), rate = d_0 + sum(lambda))
  
}

###############################################################################
# update for psi
###############################################################################

update_psi = function(Y, G, X, a, b){
  # Y: matrix with the original data, of size D \times N
  # G: mixing matrix, of size D \times K
  # X: latent variables, of size K \times N
  # a: shape parameter of prior for each psi_d^{-1}
  # b: rate parameter of prior for each psi_d^{-1}
  
  D = nrow(Y); N = ncol(Y)
  hat_E = Y - G%*%X
  
  # get matrix of parameters size K \times 2 (shape and rate)
  psi_inv_pars = matrix(NA_real_, nrow = D, ncol = 2L)
  colnames(psi_inv_pars) = c("shape", "rate")
  
  # get vector of shape parameters (length K)
  psi_inv_pars[,"shape"] = a + 0.5*N
  
  # get vector of rate parameters (length K)
  psi_inv_pars[,"rate"] = b + 0.5*rowSums(hat_E*hat_E)
  
  # draw psi_inv
  psi_inv = apply(psi_inv_pars, 1L, function(p){
    rgamma(n=1L, shape = p["shape"], rate = p["rate"])
  })
  
  # return psi
  return(1/psi_inv)
}

###############################################################################
# update for psi's parameter b
###############################################################################

update_b = function(psi, a, a_0, b_0){
  # psi: vector with the diagonal of noise variance matrix, of length D
  # a: shape parameter of prior for each psi_d^{-1}
  # a_0: shape parameter of prior for b
  # b_0: rate parameter of prior for b
  
  rgamma(n=1L, shape = a_0+a*length(psi), rate = b_0 + sum(1/psi))
  
}
