###############################################################################
# main function for mixture updates
###############################################################################

update_Z_G_X_new_lambda_new = function(alpha, lambda, Z, Y, G, X, psi, hp, d_lambda){
  # lambda: vector with precisions lambda_k for sampling g_dk, length K_+
  # Z     : (truncated) IBP matrix, of size D \times K_+.
  #         We assume that Z is fully active between columns 1 and K_+
  #         Therefore, at the end of this step we need to enforce this
  #         for the following updates to work
  # Y     : matrix with the original data, of size D \times N
  # G     : mixing matrix, of size D \times K_+
  # X     : latent variables, of size K_+ \times N
  # psi   : vector with the diagonal of noise variance matrix, of length D
  # hp    : list with hyperparameters
  # d_lambda: current rate for sampling lambda

  # get dimensions
  D = nrow(Y); K = ncol(Z)

  # useful reusable quantities
  XXT = tcrossprod(X)
  YXT = tcrossprod(Y,X)

  # loop objects
  for(d in seq_len(D)){

    # update XXT and YXT if X grew in previous iteration
    if (d > 1L && res$kappa_d > 0L){
      XXT = tcrossprod(X)
      YXT = tcrossprod(Y,X)
    }

    # update Z_d and G_d
    for(k in seq_len(K)){
      res = update_Z_dk_G_dk(d=d, k=k, lambda=lambda, Z=Z, Y=Y, G=G, X=X,
                             YXT=YXT, XXT=XXT, psi_d = psi[d])
      Z[d, k] = res$Z_dk
      G[d, k] = res$G_dk
    }

    # get new factors
    res = update_get_features_MH(d=d, Z=Z, Y=Y, G=G, X=X, lambda=lambda,
                                 alpha=alpha, psi_d=psi[d],
                                 lambda_MH=hp$lambda_MH, pi_MH=hp$pi_MH,
                                 c_lambda=hp$c, d_lambda=d_lambda)
    if(res$kappa_d > 0L){
      # grow lambda
      lambda = c(lambda, res$lambda_new)

      # grow Z
      Z_new = matrix(0L, nrow = D, ncol = res$kappa_d)
      Z_new[d, ] = rep(1L, res$kappa_d)
      Z = cbind(Z, Z_new)

      # grow G
      G_new = matrix(0, nrow = D, ncol = res$kappa_d)
      G_new[d, ] = res$g_new_T
      G = cbind(G, G_new)

      # grow X
      X = rbind(X, res$X_new)
    }
  }

  # Q: IS THIS NECESSARY? Since Z=0 for those components, and therefore
  # the corresponding G are also 0, GX doesn't change
  # A: YES!, because matrices can grow very large
  # drop inactive components them from all variables
  n_A = which(colSums(Z) == 0L)
  if(length(n_A) > 0L){
    # drop = FALSE ensures we always have matrices and not vectors
    Z = Z[, -n_A, drop = FALSE]
    G = G[, -n_A, drop = FALSE]
    X = X[-n_A, , drop = FALSE]
    lambda = lambda[-n_A]
  }

  # return
  return(list(Z=Z, G=G, X=X, lambda=lambda))

}

###############################################################################
# function to update each component Z_dk and G_dk
###############################################################################

update_Z_dk_G_dk = function(d, k, lambda, Z, Y, G, X, YXT, XXT, psi_d){
  # YXT   : tcrossprod(Y,X)
  # XXT   : tcrossprod(X)

  D = nrow(Y)

  # calculate lambda_dk = (1/psi_d) X_{k:}X_{k:}^T + lambda_k
  lambda_dk = XXT[k,k]/psi_d + lambda[k]

  # get mu_dk = 1/(lambda_dk*psi_d) \hat{E}_{d:} X_{k:}^T
  G_dk_d_XXT_k = sum(G[d, ] * XXT[,k]) - G[d,k]*XXT[k,k] # ==sum(G_star[d, ] * XXT[k,])
  hat_E_d_XT_k = YXT[d,k] - G_dk_d_XXT_k
  mu_dk = hat_E_d_XT_k / (lambda_dk*psi_d)

  # count positions in column active other than at (d,k)
  m_not_dk = sum(Z[,k]) - Z[d,k]

  if(m_not_dk == 0L){
    Z_dk = 0L
  } else{
    # get log prior odds ratio = m_{-dk} / (D - m_{-dk})
    log_prior_or_z_dk = log(m_not_dk) - log(D - m_not_dk)

    # get log lklhd ratio lk_y_z = sqrt(lambda_k/lambda_dk)exp(0.5 lambda_dk mu_dk^2)
    log_lk_y_z_dk = 0.5*(log(lambda[k]) - log(lambda_dk) + lambda_dk*(mu_dk^2))

    # get prob_z_1 = 1 - 1/(1+lk_y_z*prior_or_z)
    m = max(log_prior_or_z_dk, log_lk_y_z_dk)
    prod_lkr_prior_or = exp(m)*exp(log_prior_or_z_dk + log_lk_y_z_dk - m) # = exp(max)exp(min)
    prob_z_1 = 1 - 1/(1+prod_lkr_prior_or) # 1 - prob_z_0

    # sample Z_dk
    Z_dk = 1L*(runif(1L) < prob_z_1)
  }

  # sample G
  if(Z_dk==0L){
    G_dk = 0
  }else{
    G_dk = rnorm(n = 1L, mean = mu_dk, sd = sqrt(1/lambda_dk))
  }

  return(list(Z_dk=Z_dk, G_dk=G_dk))

}

###############################################################################
# function to get new features: affects Z, G, X, and lambda
###############################################################################

update_get_features_MH = function(d, Z, Y, G, X, lambda, alpha, psi_d, lambda_MH, pi_MH, c_lambda, d_lambda){

  # set parameters
  D = nrow(Y); N = ncol(Y)
  gamma_MH = alpha / D

  #####################################
  # get proposal \xi^*
  #####################################

  # get kappa_d proposal
  if (runif(1L) < pi_MH){
    kappa_d = 1L
  } else{
    kappa_d = rpois(n = 1L, lambda = lambda_MH*gamma_MH)
  }

  # check if kappa_d = 0 to see if we're done
  if (kappa_d == 0L){
    return(list(kappa_d=0L, lambda_new=NA_real_, g_new_T=NA_real_, X_new=NA_real_))
  }

  # get lambda proposals
  lambda_new = rgamma(n = kappa_d, shape = c_lambda, rate = d_lambda)

  # get g_new_T proposal (length kappa_d)
  g_new_T = rnorm(n = kappa_d, mean = rep(0, kappa_d), sd = sqrt(1/lambda_new))

  #####################################
  # compute acceptance ratio
  #####################################

  # get M matrix (kappa_d \times kappa_d) and associated quantities
  M = tcrossprod(g_new_T)/psi_d + diag(kappa_d)
  M_inv = solve(M)
  log_det_M = determinant(M, logarithm = TRUE)$modulus # for log(a_l)

  # get m matrix = (m_1^T,...m_n^T,...,m_N^T)^T (size N\times kappa_d)
  m = (1/psi_d) * outer(X = as.vector(Y[d,]-G[d,] %*% X), # E_d:, length N
                        Y = as.vector(M_inv %*% g_new_T)) # length kappa_d

  # get a_l efficiently, which needs trace(mMm^T)
  MmT = tcrossprod(M, m) # Mm^T (size kappa_d \times N)

  # fast trace of product of matrices
  # see https://stackoverflow.com/a/17951533/5443023
  trace_mMmT = sum(t(m) * MmT) # == sum(diag(m %*% MmT))

  # compute a_l using log for improved stability
  a_l = exp((-N/2)*log_det_M + 0.5*trace_mMmT)

  # compute a_p
  a_p_den = (1-pi_MH)*dpois(x=kappa_d, lambda=lambda_MH*gamma_MH) + pi_MH*(kappa_d==1L)
  a_p = dpois(kappa_d, gamma_MH) / a_p_den

  # compute our correction a_c
  a_c = exp(log(1-pi_MH) - gamma_MH*(lambda_MH - 1))

  # compute ratio
  # TODO: do this calculation in log space: min(0, sum of log terms).
  a_ratio = a_l*a_p*a_c

  # reject?
  if(runif(1L) > a_ratio){
    return(list(kappa_d=0L, lambda_new=NA_real_, g_new_T=NA_real_, X_new=NA_real_))
  }

  # proposal accepted!
  # need to draw X'
  U=chol(M_inv) # pre-calculate cholesky decomp of variance matrix
  X_new = apply(m, 1L, function(m_n){
    # get one sample from N(m_n, M^{-1})
    m_n + as.vector(matrix(rnorm(n = kappa_d), ncol=kappa_d) %*% U)
  }) # result is matrix of size kappa_d \times N (append below current X)

  # return
  return(list(kappa_d=kappa_d, lambda_new=lambda_new, g_new_T=g_new_T, X_new=X_new))

}

