###############################################################################
# main function for updating latent vars
###############################################################################

update_X = function(Y, G, psi){
  # Y     : matrix with the original data, of size D \times N
  # G     : mixing matrix, of size D \times K_+
  # psi   : vector with the diagonal of noise variance matrix, of length D
  
  # get common Lambda precision matrix
  K = ncol(G); N = ncol(Y)
  Lambda = crossprod(G,
                     sweep(x      = G, # dim D \times K
                           MARGIN = 1L, # each row
                           STATS  = 1/psi, # length D
                           FUN    = "*")
           )
  Lambda = Lambda + diag(K)
  Lambda_inv = solve(Lambda)
  
  # get matrix mu = (mu_1,...,mu_N) of size K \times N
  mu = tcrossprod(Lambda_inv, G) %*% sweep(x      = Y, # dim D \times N
                                           MARGIN = 1L, # each row
                                           STATS  = 1/psi, # length D
                                           FUN    = "*")
  
  # draw X
  U=chol(Lambda_inv) # pre-calculate cholesky decomp of variance matrix
  X = apply(mu, 2L, function(mu_n) {
    # get one sample from N(mu_n, Lambda^{-1})
    mu_n + as.vector(matrix(rnorm(n = K), ncol = K) %*% U)
  }) 
  
  # when K>1, result is matrix of size K \times N
  # however, it's a vector in case K=1 => need to reformat as matrix
  return(matrix(X, nrow = K, ncol = N))
  
}