#' @title poisson_likelihood
#' @description negative log-likelihood function of the three-arm study with poisson endpoints
#' @param rate_exp rate in the experimental treatment group
#' @param rate_ref rate in the reference treatment group
#' @param rate_pla rate in the placebo treatment group
#' @param n_exp sample size in the experimental treatment group
#' @param n_ref sample size in the reference treatment group
#' @param n_pla sample size in the placebo treatment group
#' @param sum_exp sum of results of the experimental treatment group
#' @param sum_ref sum of results of the reference treatment group
#' @param sum_pla sum of results of the placebo treatment group
#' @return numeric
#' @keywords internal
poisson_likelihood <- function(theta, n_exp, n_ref, n_pla, sum_exp, sum_ref, sum_pla) {
  rate_exp <- theta[1]
  rate_ref <- theta[2]
  rate_pla <- theta[3]
  out <- -(n_exp * rate_exp + n_ref * rate_ref + n_pla * rate_pla) +
    sum_exp * log(rate_exp) +
    sum_ref * log(rate_ref) +
    sum_pla * log(rate_pla)

  (-1)*out
}

#' @title grad_poisson_likelihood
#' @description gradient of negative log-likelihood function for poisson endpoints
#' @keywords internal
grad_poisson_likelihood <- function(theta, n_exp, n_ref, n_pla, sum_exp, sum_ref, sum_pla) {
  out <- numeric(3)
  out[1] <- n_exp - sum_exp / theta[1]
  out[2] <- n_ref - sum_ref / theta[2]
  out[3] <- n_pla - sum_pla / theta[3]
  out
}

#' @title negbin_likelihood
#' @description negative log-likelihood function of the three-arm study with negative binomial endpoints
#' @param theta parameter vector (rate_exp, rate_ref, rate_pla)
#' @param x_exp results of the experimental treatment group
#' @param x_ref results of the reference treatment group
#' @param x_pla results of the placebo treatment group
#' @return numeric
#' @keywords internal
negbin_likelihood <- function(theta, x_exp, x_ref, x_pla){
  rate_exp <- theta[1]
  rate_ref <- theta[2]
  rate_pla <- theta[3]
  shape <- theta[4]
  lh_exp <- sum( dnbinom(x_exp, mu = rate_exp, size = 1 / shape, log = TRUE) )
  lh_ref <- sum( dnbinom(x_ref, mu = rate_ref, size = 1 / shape, log = TRUE) )
  lh_pla <- sum( dnbinom(x_pla, mu = rate_pla, size = 1 / shape, log = TRUE) )

  (lh_exp + lh_ref + lh_pla)
}




#' @title binary_likelihood_rest
#' @description negative log-likelihood function of the three-arm study with binary endpoints
#' @param theta vector of rates
#' @param n_exp sample size in the experimental treatment group
#' @param n_ref sample size in the reference treatment group
#' @param n_pla sample size in the placebo treatment group
#' @param sum_exp sum of results of the experimental treatment group
#' @param sum_ref sum of results of the reference treatment group
#' @param sum_pla sum of results of the placebo treatment group
#' @param h Function used in hypothesis
#' @param h_inv Inverse function of h
#' @return numeric
#' @keywords internal
binary_likelihood_rest <- function(theta, Delta, n_exp, n_ref, n_pla, sum_exp, sum_ref, sum_pla, h, h_inv) {
  if (is.null(h) || is.null(h_inv)) {
    h <- identity
    h_inv <- identity
  }

  rate_ref <- theta[1]
  rate_pla <- theta[2]
  rate_exp <- h_inv(Delta * h(rate_ref) + (1-Delta) * h(rate_pla))

  if ( any(c(rate_exp, rate_ref, rate_pla) <= 0) || any(c(rate_exp, rate_ref, rate_pla) >= 1))  {
    return(99999)
  }

  out <- sum_exp * log(rate_exp) + (n_exp - sum_exp) * log(1 - rate_exp) +
    sum_ref * log(rate_ref) + (n_ref - sum_ref) * log(1 - rate_ref) +
    sum_pla * log(rate_pla) + (n_pla - sum_pla) * log(1 - rate_pla)

  (-1)*out
}
