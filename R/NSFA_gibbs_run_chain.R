###############################################################################
# initialize the chain
###############################################################################

init_chain = function(Y, hp){
  # Y  : matrix with the original data, of size D \times N
  # hp : hyperparameters

  D = nrow(Y); N = ncol(Y)

  # init alpha
  alpha = rgamma(n = 1L, shape = hp$e, rate = hp$f)

  # init Z: use finite IBP for smallest set possible -> K=2
  K = 2L
  pi_vec = rbeta(n=K, alpha/K, 1)
  Z = do.call(cbind,
              lapply(pi_vec, function(pi_k) {
                1L * (runif(n = D) < pi_k)
              }))

  # # init Z: smallest, densest matrix possible -> K=2
  # K = 2L
  # Z = matrix(1L, nrow = D, ncol = K)

  # init lambda
  d = rgamma(n=1L, shape = hp$c_0, rate = hp$d_0) # sample rate param
  lambda = rgamma(n=K, shape = hp$c, rate = d)

  # init G
  G = do.call(cbind,
              lapply(lambda, function(l_k){
                rnorm(n=D, mean = 0, sd = sqrt(1/l_k))
              }))
  G = G*Z # block inactive components

  # init X
  X = matrix(rnorm(N*K), nrow = K, ncol = N)

  # init psi
  b = rgamma(n=1L, shape = hp$a_0, rate = hp$b_0) # sample rate param
  psi = 1/(rgamma(n=D, shape = hp$a, rate = b))

  # return
  list(alpha=alpha, Z=Z, d=d, lambda=lambda, G=G, X=X, b=b, psi=psi)

}

###############################################################################
# main function for running the chain
###############################################################################

#' Gibbs Sampler for the Non-parametric Sparse Factor Analysis (NSFA) Model
#'
#' This function runs a Gibbs sampler to obtain posterior samples from the NSFA
#' model by Knowles & Ghahramani (2011). This is my personal implementation, which
#' contains some corrections with respect to the model presented in the original
#' article.
#'
#' The model was originally proposed in Knowles & Ghahramani (2011). However, you
#' can read my own report evaluating this model
#' \href{https://miguelbiron.github.io/docs/STAT548_report_3_NSFA.pdf}{here}.
#'
#' The following hyperparameters can be supplied by the user through \code{hp}.
#' Default values are in parenthesis and tend to work well.
#' \itemize{
#'   \item \code{e}: shape parameter of prior for alpha (1).
#'   \item \code{f}: rate parameter of prior for alpha (1).
#'   \item \code{a}: shape parameter of prior for psi^{-1} (1).
#'   \item \code{a_0}: shape parameter of prior for b (1).
#'   \item \code{b_0}: rate parameter of prior for b (1).
#'   \item \code{c}: shape parameter of prior for lambda (1).
#'   \item \code{c_0}: shape parameter of prior for d (1).
#'   \item \code{d_0}: rate parameter of prior for d (1).
#'   \item \code{lambda_MH}: lambda parameter for kappa_d proposal (10).
#'   \item \code{pi_MH}: pi parameter for kappa_d proposal (0.1).
#' }
#'
#' If chain keeps restarting, try increasing \code{lambda_MH} and/or \code{pi_MH}
#' to get more aggressive proposals for growing Z.
#'
#' @param S number of iterations for the chain.
#' @param Y matrix with the observed data, of size \code{D x N}.
#' @param hp list of hyperparameters. See Details for a description and default
#' values used.
#' @param Y_test matrix with test set for assessing out-of-sample performance.
#' Number of columns must be same as \code{Y} (\code{N}).
#' @param R function returns last \code{R} values in the chain. Default is 1.
#' @param verbose integer. Function prints every \code{verbose} iterations. The
#' default value of 0 prints no output.
#'
#' @return List with two lists:
#' \describe{
#'   \item{chain}{list with last \code{R} samples of the chain.}
#'   \item{stats}{list of traces of train set log-likelihood (ll), test set ll, and K.}
#' }
#'
#' @references
#' Knowles, D., & Ghahramani, Z. (2011). Nonparametric Bayesian sparse factor
#' models with application to gene expression modeling.
#' \emph{The Annals of Applied Statistics}, 1534-1552.
#'
#' @export
NSFA_gibbs = function(S, Y, hp = NULL, Y_test = NULL, R = 1L, verbose = 0L){

  stopifnot(is.matrix(Y) && R <= S) # sanity checks

  # check hp
  if(is.null(hp)){
    hp = vector(mode = "list", length = 10L)
    names(hp) = c("e", "f", "a", "a_0", "b_0", "c", "c_0", "d_0", "lambda_MH", "pi_MH")

    # all hyperparameters equal to 1
    # from: https://github.com/davidaknowles/nsfa/code/defaultsettings.m
    hp$e = hp$a = hp$a_0 = hp$c = hp$c_0 = hp$f = hp$b_0 = hp$d_0 = 1

    # MH proposal. lambda from same source, pi from paper.
    hp$lambda_MH = 10; hp$pi_MH = 0.1
  } else{
    stopifnot(length(hp) == 10L) # sanity check
  }

  # initialize chain
  chain = init_chain(Y, hp)
  return_chain = vector(mode = "list", length = R)

  # storage of useful statistics
  stats = vector(mode = "list")
  stats$ll_Y = numeric(S) # storage for log likelihood
  stats$ll_Y[1L] = loglik_Y(Y = Y,G = chain$G,X = chain$X,psi = chain$psi)
  stats$K_vec = integer(S)
  stats$K_vec[1L] = ncol(chain$Z)
  if(!is.null(Y_test)){
    stopifnot(ncol(Y_test) == ncol(Y)) # sanity check
    D_test = nrow(Y_test)
    stats$ll_Y_test = numeric(S) # storage for test set log likelihood
    ho_pars = held_out_eval(chain = chain, hp = hp, D_test = D_test)
    stats$ll_Y_test[1L] = loglik_Y(Y=Y_test,G=ho_pars$G_test, X=chain$X, psi=ho_pars$psi_test)
  }

  # calculate reusable quantities
  D = nrow(Y)
  H_D = sum(1/(seq_len(D))) # D-th harmonic number

  # loop
  tryCatch({
    for(s in 2L:S){

      # print status
      if(verbose > 0L && (s %% verbose == 0L)){
        cat(sprintf("Iteration %d: pre-K=%d", s, stats$K_vec[s-1L]))
      }

      # update mixture components
      res = update_Z_G_X_new_lambda_new(
        alpha    = chain$alpha,
        lambda   = chain$lambda,
        Z        = chain$Z,
        Y        = Y,
        G        = chain$G,
        X        = chain$X,
        psi      = chain$psi,
        hp       = hp,
        d_lambda = chain$d
      )
      if(length(res$lambda) > 0L){
        chain$lambda = res$lambda
        chain$Z = res$Z
        chain$G = res$G
        chain$X = res$X # only adds new rows if kappa_d>0 was accepted
      } else{
        # everything was erased, need to restart
        chain = init_chain(Y=Y, hp=hp)
        if(verbose > 0L && (s %% verbose == 0L)){
          cat(", chain restarted!\n")
        }
        next()
      }

      # update X
      chain$X = update_X(
        Y   = Y,
        G   = chain$G,
        psi = chain$psi
      )

      # update alpha
      K = ncol(chain$Z) # get current K
      chain$alpha = update_alpha(K=K, H_D=H_D, e=hp$e, f=hp$f)

      # update lambda
      chain$lambda = update_lambda(
        G = chain$G,
        Z = chain$Z,
        c = hp$c,
        d = chain$d
      )

      # update rate for lambda d
      chain$d = update_d(
        lambda = chain$lambda,
        c      = hp$c,
        c_0    = hp$c_0,
        d_0    = hp$d_0
      )

      # update psi
      chain$psi = update_psi(
        Y = Y,
        G = chain$G,
        X = chain$X,
        a = hp$a,
        b = chain$b
      )

      # update rate for sampling psi (b)
      chain$b = update_b(
        psi = chain$psi,
        a   = hp$a,
        a_0 = hp$a_0,
        b_0 = hp$b_0
      )

      # store likelihood and K
      stats$ll_Y[s] = loglik_Y(Y = Y,G = chain$G,X = chain$X,psi = chain$psi)
      stats$K_vec[s] = ncol(chain$Z)
      if(!is.null(Y_test)){
        ho_pars = held_out_eval(chain = chain, hp = hp, D_test = D_test)
        stats$ll_Y_test[s] = loglik_Y(Y=Y_test,G=ho_pars$G_test, X=chain$X, psi=ho_pars$psi_test)
      }

      # collect current point of chain in return_chain
      if(R > 1L && s > (S-R)){
        return_chain[[s-(S-R)]] = chain
      }

      # print status
      if(verbose > 0L && (s %% verbose == 0L)){
        cat(sprintf(", post-K=%d, loglik=%.2f", stats$K_vec[s], stats$ll_Y[s]))
        if(!is.null(Y_test)){
          cat(sprintf(", test_loglik=%.2f.\n", stats$ll_Y_test[s]))
        }else{
          cat(".\n")
        }
      }
    }

    # set iterations names for return_chain
    names(return_chain) = seq.int(S-R+1L, S, by=1L)

  }, error = function(e){
    cat("\n\n")
    print(e)
    warning("Error found. Inspect results to understand why, or increase verbosity.")
  }, finally = {return(list(chain=return_chain, stats=stats))})

}


