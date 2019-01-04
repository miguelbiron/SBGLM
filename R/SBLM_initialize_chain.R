initialize_from_prior = function(hp, P){
  
  # sample sigma_2_b and sigma_2_e
  # Ideally, we want sigma_2_e << sigma_2_b (high SNR ratio)
  sigma_2_b = 1/rgamma(n=1L, shape=hp$a_0_b, rate=hp$b_0_b)
  sigma_2_e = 1/rgamma(n=1L, shape=hp$a_0_e, rate=hp$b_0_e)
  
  # sample pi
  pi_z = rbeta(n=1L, shape1 = hp$s_1_pi, shape2 = hp$s_2_pi)
  
  # sample switches
  z = sample(c(TRUE, FALSE), size=P, replace=TRUE, prob=c(pi_z, 1-pi_z))
  A = which(z == TRUE) # active set
  
  # sample beta
  beta = rnorm(n = P, mean = 0, sd = sqrt(sigma_2_b))
  
  # return
  list(
    A         = A,
    sigma_2_b = sigma_2_b,
    sigma_2_e = sigma_2_e,
    beta      = beta,
    pi_z      = pi_z
  )
  
}

initialize_from_OLS = function(y, X){
  
  fit_ols = lsfit(x=X, y=y, intercept = FALSE)
  beta = fit_ols$coefficients
  sigma_2_e = mean(fit_ols$residuals^2)
  A = seq_len(ncol(X))
  sigma_2_b = mean(fit_ols$coefficients^2)
  
  list(
    A         = A,
    sigma_2_b = sigma_2_b,
    sigma_2_e = sigma_2_e,
    beta      = beta,
    pi_z      = mean(abs(beta)>0.000001)
  )
}