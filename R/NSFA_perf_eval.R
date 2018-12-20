###############################################################################
# log-likelihood of the data (up to additive constant)
###############################################################################

loglik_Y = function(Y, G, X, psi){
  N = ncol(Y)
  E = Y - G %*% X
  -0.5*(N*sum(log(psi)) + sum(rowSums(E*E) / psi))
}

###############################################################################
# out-of-sample evaluation
# sample psi, Z and G for held-out data
###############################################################################

held_out_eval = function(chain, hp, D_test){
  # D_test = nrow(Y_test)
  
  D = nrow(chain$Z) # train set D
  K = nrow(chain$X)
  
  #####################################
  # sample Z_test
  #####################################
  
  # given K, we can use the finite version of the IBP to sample Z
  # need to sample pi_1...pi_K | alpha, K, Z_train
  m = colSums(chain$Z) # vector containing the m_k's
  shape = matrix(NA_real_, nrow = K, ncol = 2L)
  shape[,1L] = chain$alpha/K + m # vector containing the shape1 parameters (length K)
  shape[,2L] = D - m + 1 # vector containing the shape2 parameters (length K)
  
  # sample pi (vector length K)
  pi = apply(shape, 1L, function(s){
    rbeta(n=1L, shape1 = s[1], shape2 = s[2])
  })
  
  # sample Z_test
  Z_test = do.call(cbind,
                   lapply(pi, function(pi_k){
                     1L*(runif(n = D_test) < pi_k)
                   }))
  
  #####################################
  # sample G and psi
  #####################################
  
  # sample G
  G_test = do.call(cbind,
              lapply(chain$lambda, function(l_k){
                rnorm(n=D_test, mean = 0, sd = sqrt(1/l_k))
              }))
  G_test = G_test*Z_test # block inactive components
  
  # sample psi
  psi_test = 1/(rgamma(n=D_test, shape = hp$a, rate = chain$b))
  
  return(list(Z_test=Z_test, G_test=G_test, psi_test=psi_test))
  
}

# # test
# ho_pars = held_out_eval(chain, hp, nrow(Y))
# loglik_Y(Y = Y,G = ho_pars$G_test, X = chain$X, psi = ho_pars$psi_test)
