#' @title Power related calculations for three-arm clinical trials
#' @description Compute power, sample size, or level of significance for
#' Wald-type test for non-inferiority or superiority of the experimental treatment
#' versus reference treatment with respect to placebo.
#' @details If the individual group sample sizes, i.e. \code{n*allocation}
#' are not natural number, the parameters \emph{n} and \emph{allocation}
#' will be re-calculated.
#' @param experiment a numeric vector specifying the parameters of the experimental
#' treatment group in the alternative hypothesis
#' @param reference a numeric vector specifying the parameters of the reference
#' treatment group in the alternative hypothesis
#' @param placebo a numeric vector specifying the parameters of the placebo
#' treatment group in the alternative hypothesis
#' @param Delta a numeric value specifying the non-inferiority/superiority margin
#' @param sig_level A numeric value specifying the significance level (type I error probability)
#' @param power A numeric value specifying the target power (1 - type II error probability)
#' @param n The total sample size. Needs to be at least 7.
#' @param allocation A (non-empty) vector specifying the sample size allocation (nExp/n, nRef/n, nPla/n)
#' @param distribution A character specifying the distribution of the endpoints. Must
#'        must be either of \code{"binary"}, \code{"poisson"}, \code{"negbin"}, \code{"exponential"}, \code{"normal"}
#' @param ... Further arguments. See details.
#' @details
#' The additional parameter \code{var_estimation} is a character string
#' specifying how the variance for the Wald-type test statistic is estimated
#' in the Poisson and negative binomial model. Must be \emph{RML} for restricted
#' maximum-likelihood, or \emph{ML} for unrestricted maximum-likelihood
#' @return A list with class "power.htest" containing the following components:
#' \item{n}{The total sample size}
#' \item{power}{A numeric value specifying the target power}
#' \item{Delta}{A numeric value specifying the non-inferiority or superiority margin. }
#' \item{sig.level}{A character string specifying the significance level}
#' \item{type}{A character string indicating what type of Wald-type test will be performed}
#' \item{allocation}{A vector with the sample size allocation (nExp/n, nRef/n, nPla/n)}
#' \item{sig.level}{The significance level (Type I error probability)}
#' \item{nExp}{A numeric value specifying the number of sample in the experimental treatment group}
#' \item{nRef}{A numeric value specifying the number of sample in the reference treatment group}
#' \item{nPla}{A numeric value specifying the number of sample in the placebo treatment group}
#' @examples
#' power_RET(experiment = 15, reference = 17, placebo = 20,
#'          Delta = 0.8, sig_level = 0.025, power = 0.8,
#'          allocation = c(1, 1, 1) / 3,
#'          var_estimation = "RML",
#'          distribution = "poisson")
#' @export
#' @keywords power waldtype samplesize
power_RET <- function(experiment, reference, placebo,
                      Delta,
                      sig_level = NULL,
                      power = NULL,
                      n = NULL,
                      allocation = c(1/3, 1/3, 1/3),
                      distribution = NULL,
                      ...){

  arguments <- names(formals())
  arguments_dots <- list(...)

  check_missing(args = arguments[!(arguments %in% c("n", "sig_level", "power", "distribution", "..."))])
  check_RET_arguments(sig.level = sig_level, power = power,
                      Delta = Delta, n = n, allocation = allocation)

  if (is.null(distribution)) {
    warning("Input parameters are considered as group means and variances", call. = FALSE)
  } else{
    match.arg(arg = distribution,
              choices = c("poisson", "negbin", "normal", "exponential", "binary"),
              several.ok = FALSE)
  }

  if( sum(sapply(list(n, power, sig_level), is.null)) !=  1 )
    stop("Exactly one of 'n', 'power', and 'sig_level' must be NULL.")

  x <- init_power(para_exp = experiment,
                  para_ref = reference,
                  para_pla = placebo,
                  distribution = distribution)

  # Adjust sample size and allocation if necessary
  w <- allocation
  if (!is.null(n)) {
    n_exp <- round(n * w[1])
    n_ref <- round(n * w[2])
    n_pla <- round(n * w[3])
    n_group <- c(n_exp, n_ref, n_pla)
    if (!identical(n_group / sum(n_group), w)) {
      n <-  sum(n_group)
      w <- n_group / n
      warning("Allocation not suitable for current sample size.  'allocation' and 'n' will be adjusted.")
    }
    else if (sum(n_group) != n) {
      n <-  sum(n_group)
      warning("'n' not suitable for the defined allocation and will be adjusted.")
    }
  }

  calc_power_ret(x = x,
                 Delta = Delta,
                 allocation = w,
                 n = n,
                 power = power,
                 sig_level = sig_level,
                 distribution = distribution, ...)

}


