#' @title Wald-type test for three-arm trials
#' @description Wald-type test for superiority/non-inferiority of the
#' experimental treatment versus reference treatment with respect to placebo.
#' @param xExp A (non-empty) numeric vector of data values
#' from the experimental treatment group.
#' @param xRef A (non-empty) numeric vector of data values
#' from the reference treatment group.
#' @param xPla A (non-empty) numeric vector of data values
#' from the placebo group.
#' @param Delta A numeric value specifying the non-inferiority or superiority margin.
#' Is between 0 and 1 in case of non-inferiority and larger than 1 in case of superiority.
#' @param ... Other named arguments such as \code{distribution}, \code{var_estimation}.
#' See details for more information.
#' @details
#' Additional parameters include \code{distribution} and \code{var_estimation}. \cr
#' The parameter \code{distribution} is a character string and indicates
#' whether a parametric model should be used. If not specified retention of
#' effect hypothesis is tested using sample means and variances.
#' The following options exist:
#' \code{"poisson"} (Poisson distribution),
#' \code{"negbin"} (negative binomial distribution),
#' \code{"normal"} (normal distribution),
#' \code{"exponential"} (censored exponential).
#' \code{"nonparametric"} (non-parametric).
#' If the parameter \code{distribution} is not specified
#' the effect and the variance for the test statistic are estimated
#' by the sample means and sample variances.\cr
#' The parameter \code{var_estimation} defines how the variance is estimated
#' in the parametric models \code{"poisson"} and \code{"negbin"}.
#' The following options exist:
#' \code{RML} for the restricted maximum-likelihood estimator
#' and \code{ML} (default) for the unrestricted maximum-likelihood estimator.
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{The value of the Wald-type test statistic.}
#' \item{p.value}{The p-value for the Wald-type test.}
#' \item{method}{A character string indicating what type of Wald-type-test was performed.}
#' \item{estimate}{The estimated rates for each of the group as well as the maximum-likelihood estimator for the shape parameter.}
#' \item{sample.size}{The total number of data points used for the Wald-type test.}
#' @examples
#' # Negative binomially distributed endpoints
#' # Test for non-inferiority test. lambda_P=8, lambda_R = 4, lambda_E = 5, and phi = 1
#' # Delta = (lambda_P-lambda_E)/(lambda_P-lambda_R)
#' xExp <- rnbinom(60, mu = 5, size = 1)
#' xRef <- rnbinom(40, mu = 4, size = 1)
#' xPla <- rnbinom(40, mu = 8, size = 1)
#' Delta <- (8-5) / (8-4)
#' test_RET(xExp, xRef, xPla, Delta, var_estimation = 'RML', distribution = "negbin")
#' test_RET(xExp, xRef, xPla, Delta, var_estimation = 'ML', distribution = "negbin")
#'
#' # Poisson distributed endpoints
#' # Test for non-inferiority test. lambda_P=8, lambda_R = 4, lambda_E = 5
#' # Delta = (lambda_P-lambda_E)/(lambda_P-lambda_R)
#' xExp <- rpois(60, lambda = 5)
#' xRef <- rpois(40, lambda = 4)
#' xPla <- rpois(40, lambda = 8)
#' Delta <- (8-5) / (8-4)
#' test_RET(xExp, xRef, xPla, Delta, var_estimation = 'RML', distribution = "poisson")
#' test_RET(xExp, xRef, xPla, Delta, var_estimation = 'ML', distribution = "poisson")
#'
#' # Censored exponential distributed endpoints
#' # Test for non-inferiority test. lambda_P=3, lambda_R = 1, lambda_E = 2
#' # Probability for uncensored observation: 0.9
#' # Delta = (lambda_P-lambda_E)/(lambda_P-lambda_R)
#' x_exp <- matrix(c(rexp(40, rate = 1/2), rbinom(40, size = 1, prob = 0.9)),
#'                  ncol = 2, byrow = FALSE)
#' x_ref <- matrix(c(rexp(40, rate = 1/1), rbinom(40, size = 1, prob = 0.9)),
#'                  ncol = 2, byrow = FALSE)
#' x_pla <- matrix(c(rexp(40, rate = 1/3), rbinom(40, size = 1, prob = 0.9)),
#'                  ncol = 2, byrow = FALSE)
#' Delta <- log(2/3) / log(1/3)
#' test_RET(xExp = x_exp,
#'                  xRef = x_ref,
#'                  xPla = x_pla,
#'                  Delta = Delta,
#'                  distribution = "exponential")
#' @references
#' I. Pigeot, J. Schaefer, J. Roehmel, D. Hauschke. (2008).
#' \emph{Assessing non-inferiority of a new treatment in a three-arm clinical trial including a placebo.}
#' Statistics in Medicine.  30(6):883-99.
#'
#' M. Hasler, R. Vonk, and LA. Hothorn. (2008).
#' \emph{Assessing non-inferiority of a new treatment in a three-arm trial in the presence of heteroscedasticity.}
#' Statistics in Medicine, 27(4):490-503.
#'
#' M. Mielke and A. Munk. (2009).
#' \emph{The assessment and planning of non-inferiority trials for retention of effect hypotheses-towards a general approach.}
#' arXiv preprint arXiv:0912.4169.
#'
#' T. Muetze, A. Munk, and T. Friede. (2016).
#' \emph{Design and analysis of three-arm trials with negative binomially distributed endpoints.}
#' Statistics in Medicine, 35(4):505-521.
#' @seealso \code{\link{power_RET}}
#' @export
#' @keywords test waldtype
test_RET <- function (xExp, xRef, xPla, Delta, ...) {

  if (missing(Delta) || is.null(Delta) || (Delta <= 0) || (!is.numeric(Delta))) {
    stop("'Delta' must be numeric and positive.")
  }

  arguments <- list(...)
  if (is.null(arguments$distribution)) {
    warning("Effect and test statistic variances are estimated using the mean and sample variances.", call. = FALSE)
  } else{
    match.arg(arg = arguments$distribution,
              choices = c("poisson", "negbin", "normal", "exponential", "nonparametric", "binary"),
              several.ok = FALSE)
  }

  data_name <- paste(c(deparse(substitute(xExp)), ', ',
                       deparse(substitute(xRef)), ', and ',
                       deparse(substitute(xPla)),
                       collapse = ''))

  x <- init_test(data_exp = xExp,
                 data_ref = xRef,
                 data_pla = xPla,
                 distribution = arguments$distribution)

  estimate_out <- estimate(x = x, Delta = Delta, ...)
  x$effect <- estimate_out$effect
  x$teststat_var <- estimate_out$teststat_var
  x$group_response <- estimate_out$group_response

  calc_test_ret(x = x, Delta = Delta, data_name = data_name, ...)
}
