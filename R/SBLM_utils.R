# get Sigma and sqrt(Sigma) from Sigma^{-1}
get_sigma_sqrt_sigma = function(A){
  # get eigen decomposition of A
  eigen_A = eigen(x = A, symmetric = TRUE)
  
  # recover Sigma = A^{-1} = Q %*% D^{-1} %*% Q^T
  # this is a fast version of the above
  Sigma = tcrossprod(
    eigen_A$vectors,
    sweep(
      x      = eigen_A$vectors,
      MARGIN = 2L,
      STATS  = 1 / eigen_A$values,
      FUN    = "*"
    )
  )
  
  # get Sigma^{1/2} = Q %*% D^{-1/2} Q^T
  sqrt_Sigma = tcrossprod(eigen_A$vectors,
                          sweep(
                            x      = eigen_A$vectors,
                            MARGIN = 2L,
                            STATS  = 1 / sqrt(eigen_A$values),
                            FUN    = "*"
                          ))
  
  return(list(Sigma=Sigma, sqrt_Sigma=sqrt_Sigma))
  
}

# sample multivariate normal given sqrt(Sigma)
rmvnorm = function(n, mu, sqrt_Sigma){
  stopifnot(length(mu) == ncol(sqrt_Sigma))
  P = ncol(sqrt_Sigma)
  if(n>1L){
    sweep(x      = matrix(rnorm(n*P), ncol = P) %*% sqrt_Sigma,
          MARGIN = 1L,
          STATS  = mu,
          FUN    = "+")
  } else{
    mu + sqrt_Sigma %*% rnorm(P)
  }
  
}
