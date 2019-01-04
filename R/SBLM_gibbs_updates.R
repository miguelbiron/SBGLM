update_sigma_2_e = function(beta, y, X, A, a_0_e, b_0_e){
  
  # update a_n
  N = length(y)
  a_n = a_0_e + 0.5*N
  
  # get predicted mean
  if(length(A) > 1){
    X_beta = as.vector(X[,A] %*% beta[A])
  } else if(length(A) == 1){
    X_beta = X[,A]*beta[A]
  } else{
    X_beta = rep(0, N)
  }
  
  # update b_n
  b_n = b_0_e + 0.5*sum((y - X_beta)^2)
  
  # sample sigma_2_e
  return(1/rgamma(n=1L, shape=a_n, rate=b_n))
}

update_sigma_2_b = function(beta, a_0_b, b_0_b){
  
  # update b_n
  P = length(beta)
  a_n = a_0_b + 0.5*P
  b_n = b_0_b + 0.5*sum(beta^2)
  
  # sample sigma_2_b
  return(1/rgamma(n=1L, shape=a_n, rate=b_n))
}

update_pi = function(A, P, s_1_pi, s_2_pi){
  rbeta(n=1L, shape1 = s_1_pi + length(A), shape2 = s_2_pi + P - length(A))
}

update_beta = function(sigma_2_b, sigma_2_e, y, X, A){
  P = ncol(X)
  beta = numeric(P)
  A_not = setdiff(seq_len(P), A)
  
  # sample beta_A_not
  beta[A_not] = rnorm(n=length(A_not), sd=sqrt(sigma_2_b))
  
  # sample beta_A
  if(length(A) > 0){
    Sigma_inv = crossprod(X[, A]) / sigma_2_e + diag(length(A)) / sigma_2_b
    res = get_sigma_sqrt_sigma(Sigma_inv)
    mu_n = as.vector(res$Sigma %*% crossprod(X[, A], y)) / sigma_2_e
    beta[A] = rmvnorm(n = 1L, mu = mu_n, sqrt_Sigma = res$sqrt_Sigma)
  }
  
  return(beta)
}

update_switches = function(beta, A, sigma_2_e, pi_z, XTy, XTX){
  # XTX = crossprod(X)
  # XTy = as.vector(crossprod(X,y))
  
  P = ncol(XTX)
  z = logical(P)
  z[A] = TRUE
  
  # precompute useful quantities
  z_0 = z
  pre_log_r = log(pi_z)-log(1-pi_z)-0.5*(1/sigma_2_e)*(beta*(-2*XTy+beta*diag(XTX)))
  
  # update switches
  for(k in seq_len(P)){
    z_0[k] = FALSE
    log_r = pre_log_r[k] - (beta[k]/sigma_2_e) * sum(XTX[,k] * beta * z_0)
    p_out = 1/(1 + exp(log_r))
    z[k] = (runif(1L) > p_out)
    z_0[k] = z[k]
  }
  
  # return set of active switches (A)
  return(which(z == TRUE))
  
}

