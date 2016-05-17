#' @title Optimal sample size for three-arm trials when analyzed with a Wald-type test
#' @description Calculate optimal sample size allocation for Wald-type test for
#' superiority or non-inferiority of the experimental treatment versus
#' reference treatment with respect to placebo
#' @param experiment a numeric vector specifying the parameters of the experimental
#' treatment group in the alternative hypothesis
#' @param reference a numeric vector specifying the parameters of the reference
#' treatment group in the alternative hypothesis
#' @param placebo a numeric vector specifying the parameters of the placebo
#' treatment group in the alternative hypothesis
#' @param Delta a numeric value specifying the non-inferiority/superiority margin
#' @param distribution a character specifying the distribution of the endpoints. Must
#'        must be either of \code{"poisson"}, \code{"negbin"}, \code{"exponential"}, \code{"normal"}
#' @param h Function measuring the efficacy; used to defined hypothesis
#' @details
#' The arguments \code{experiment}, \code{reference}, and \code{placebo} define
#' the parameters of the endpoint distribution for the respective groups: \cr
#' \code{distribution = "poisson"}: \code{experiment}, \code{reference}, and
#' \code{placebo} must have length one and define the means. \cr
#' \code{distribution = "negbin"}: \code{experiment}, \code{reference}, and
#' \code{placebo} must have length two and define the mean in the first entry
#' and the shape parameter in the second entry. \cr
#' \code{distribution = "exponential"}: \code{experiment}, \code{reference}, and
#' \code{placebo} must have length two and define the mean in the first entry
#' and the probability for an uncensored observation in the second entry. \cr
#' \code{distribution = "normal"}: \code{experiment}, \code{reference}, and
#' \code{placebo} must have length two and define the mean in the first entry
#' and the variance in the second entry.
#' @return Vector with optimal sample size allocation in the order
#' (experiment, reference, placebo)
#' @examples
#' opt_alloc_RET(experiment = 1,
#'               reference = 1,
#'               placebo = 3,
#'               Delta = 0.8,
#'               distribution = "poisson")
#'
#' @export
#' @keywords allocation waldtype
opt_alloc_RET <- function(experiment, reference, placebo, Delta, distribution, h = NULL) {

  match.arg(arg = distribution, choices = c("binary", "poisson", "negbin", "exponential", "normal"), several.ok = FALSE)

  if (is.null(h)) {
    h <- identity
  }

  if (Delta <= 0) {
    stop("Margin must be positive.")
  }
  if (any(c(missing(experiment), missing(reference), missing(placebo)))) {
    stop("None of the parameter vectors must be missing.")
  }

  length_vec <- c(length(experiment), length(reference), length(placebo))
  paras <- c(experiment, reference, placebo)


  if (distribution == "poisson") {
    if (any(length_vec != 1)) {
      stop("Only one parameter must be defined for optimal allocation calculations for poisson endpoints.")
    }
    var_exp <- experiment
    var_ref <- reference
    var_pla <- placebo
    if (any(paras <= 0)) {
      stop("Rates must be positive.")
    }
  } else if (distribution == "negbin") {
    if (any(length_vec != 2)) {
      stop("Two parameters must be defined for optimal allocation calculations for negative binomial endpoints.")
    }
    if (any(paras <= 0)) {
      stop("Rates and shape parameters must be positive.")
    }
    var_exp <- experiment[1] * (1 + experiment[1] * experiment[2])
    var_ref <- reference[1] * (1 + reference[1] * reference[2])
    var_pla <- placebo[1] * (1 + placebo[1] * placebo[2])
  } else if (distribution == "exponential") {
    if (any(length_vec != 2)) {
      stop("Two parameters must be defined for optimal allocation calculations for censored exponential endpoints.")
    }
    if (any(paras <= 0)) {
      stop("Rates and uncensore-probabilities must be positive.")
    }
    var_exp <- 1/experiment[2]
    var_ref <- 1/reference[2]
    var_pla <- 1/placebo[2]
  } else if (distribution == "normal") {
    if (any(length_vec != 2)) {
      stop("Two parameters must be defined for optimal allocation calculations for normally distributed endpoints.")
    }
    var_exp <- experiment[2]
    var_ref <- reference[2]
    var_pla <- placebo[2]
    if (any(c(var_exp, var_ref, var_pla) <= 0)) {
      stop("Variances must be positive.")
    }
  } else if (distribution == "binary") {
    if (any(length_vec != 1)) {
      stop("Parameters must have length one for optimal allocation calculations for binary endpoints.")
    }
    group_var <- numDeriv::grad(func = h, x = paras)^2 * paras * (1 - paras)
    var_exp <- group_var[1]
    var_ref <- group_var[2]
    var_pla <- group_var[3]
    if (any(c(var_exp, var_ref, var_pla) <= 0)) {
      stop("Variances must be positive.")
    }
  }

  # Optimal Allocation
  w <- c(1, Delta * sqrt(var_ref/var_exp), abs(1-Delta) * sqrt(var_pla/var_exp))
  w <- w / sum(w)
  w
}
