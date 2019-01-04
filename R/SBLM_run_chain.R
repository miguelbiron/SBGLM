###############################################################################
# main function for running the chain
###############################################################################

#' Gibbs sampler for a sparse Bayesian linear regression model
#'
#' This function runs a Gibbs sampler for a Bayesian linear regression model
#' that explicitly allows for sparse solutions, in the spirit of the spike and
#' slab prior (Mitchell and Beauchamp 1988). The difference with this model
#' is that, instead of putting a prior on the coefficients with a point mass at
#' 0, we use binary auxiliary variables to mask the effect of a variable from
#' the predicted values. More details about the model can be found in
#' \href{https://miguelbiron.github.io/}{UPDATE THIS LINK WHEN POST IS UP}.
#'
#' The following hyperparameters can be supplied by the user through \code{hp}.
#' Default values for all are 1, and tend to work well.
#' \itemize{
#'   \item \code{a_0_b}: shape parameter of Inv-Gamma prior for sigma_2_b.
#'   \item \code{b_0_b}: scale parameter of Inv-Gamma prior for sigma_2_b.
#'   \item \code{a_0_e}: shape parameter of Inv-Gamma prior for sigma_2_e.
#'   \item \code{b_0_e}: scale parameter of Inv-Gamma prior for sigma_2_e.
#'   \item \code{s_1_pi}: shape 1 parameter for Beta prior of pi_z.
#'   \item \code{s_2_pi}: shape 2 parameter for Beta prior of pi_z.
#' }
#'
#' Note that the default values for the prior of pi_z imply a uniform distribution
#' on (0,1). If the number of columns in \code{X} is large, and you suspect that
#' beta should be highly sparse, then reflecting this in the prior of pi_z would
#' speed up computation considerably.
#'
#'
#' @param X design matrix, of size \code{N x P}.
#' @param y response vector of length \code{N}.
#' @param hp list of hyperparameters. See Details for a description and default
#' values used.
#' @param S number of iterations for the chain.
#' @param verbose integer. Function prints every \code{verbose} iterations. The
#' default value of 0 prints no output.
#'
#' @return A nested list with \code{S} elements, each one being a list of the
#' values of the latent variables at every iteration.
#'
#' @references
#' Mitchell, Toby J, and John J Beauchamp. 1988. Bayesian Variable Selection in
#' Linear Regression. \emph{Journal of the American Statistical Association,
#' 83} (404). Taylor & Francis Group: 1023–32.
#'
#' @export
sblm_gibbs = function(X, y, hp = NULL, S, verbose = 0L){

  stopifnot(nrow(X)==length(y)) # sanity check

  # check if hp is provided
  if(is.null(hp)){
    hp = list()
    hp$a_0_b = hp$a_0_e = hp$b_0_b = hp$b_0_e = hp$s_1_pi = hp$s_2_pi = 1
  } else{
    stopifnot(length(hp) == 6L) # sanity check
  }

  # set problem parameters
  N = length(y)
  P = ncol(X)

  # pre-calculate required quantities
  XTX = crossprod(X)
  XTy = as.vector(crossprod(X,y))

  # initialize chain
  init = initialize_from_prior(hp=hp, P=P)

  # define storage for values of the chain
  chain = replicate(n = S, expr = init, simplify = FALSE)

  # run chain
  cat("Chain initialized. Beginning iterations!\n\n")
  tryCatch({
    for (s in 2L:S) {

      # print information
      if(verbose > 0L && (s %% verbose == 0L)) {
        cat(
          sprintf(
            "Iteration %d: |A| = %d, sigma_2_e = %.2f, sigma_2_b = %.2f\n",
            s,
            length(chain[[s - 1L]]$A),
            chain[[s - 1L]]$sigma_2_e,
            chain[[s - 1L]]$sigma_2_b
          )
        )
      }

      # update sigma_2_b
      chain[[s]]$sigma_2_b = update_sigma_2_b(
        beta  = chain[[s - 1L]]$beta,
        a_0_b = hp$a_0_b,
        b_0_b = hp$b_0_b
      )

      # update sigma_2_e
      chain[[s]]$sigma_2_e = update_sigma_2_e(
        beta = chain[[s - 1L]]$beta,
        y    = y,
        X    = X,
        A    = chain[[s - 1L]]$A,
        a_0_e = hp$a_0_e,
        b_0_e = hp$b_0_e
      )

      # update beta
      chain[[s]]$beta = update_beta(
        sigma_2_b = chain[[s]]$sigma_2_b,
        sigma_2_e = chain[[s]]$sigma_2_e,
        y         = y,
        X         = X,
        A         = chain[[s - 1L]]$A
      )

      # update switches
      chain[[s]]$A = update_switches(
        beta      = chain[[s]]$beta,
        A         = chain[[s - 1L]]$A,
        sigma_2_e = chain[[s]]$sigma_2_e,
        pi_z      = chain[[s - 1L]]$pi_z,
        XTy       = XTy,
        XTX       = XTX
      )

      # update pi_z
      chain[[s]]$pi_z = update_pi(
        A      = chain[[s]]$A,
        P      = P,
        s_1_pi = hp$s_1_pi,
        s_2_pi = hp$s_2_pi
      )
    }
  }, warning = function(w) {
    cat("\n\n")
    print(warnings())
  }, error = function(e) {
    cat("\n\n")
    print(e)
    warning("Error found. Inspect results to understand why, or increase verbosity.")
  }, finally = {
    return(chain)
  })

}
