#' @title estimate.exponential
#' @description Parameter estimation for Wald-type test with normally distributed endpoints
#' @keywords internal
estimate.exponential <- function(x, Delta, ...) {

  h <- "log"

  rate_exp <- sum(x$data_exp[, 1]) / sum(x$data_exp[, 2])
  rate_ref <- sum(x$data_ref[, 1]) / sum(x$data_ref[, 2])
  rate_pla <- sum(x$data_pla[, 1]) / sum(x$data_pla[, 2])

  effect <- do.call(what = h, args = list(x = rate_exp)) -
    Delta * do.call(what = h, args = list(x = rate_ref)) +
    (Delta - 1)  * do.call(what = h, args = list(x = rate_pla))

  teststat_var <- x$n * (1 / sum(x$data_exp[, 2]) +
                           Delta^2 / sum(x$data_ref[, 2]) +
                           (1-Delta)^2 / sum(x$data_pla[, 2]))

  list(effect = effect,
       group_response = c(rate_exp, rate_ref, rate_pla),
       teststat_var = teststat_var)
}



#' @title calc_power_ret.exponential
#' @description Power related calculations for Wald-type test with exponentially distributed endpoints
#' @keywords internal
calc_power_ret.exponential <- function(x, Delta, allocation, n = NULL, power = NULL, sig_level = NULL, ...) {

  h <- log

  if (any(c(length(x$para_exp), length(x$para_ref), length(x$para_pla)) != 2)) {
    stop("Two parameters must be defined for power related calculations for exponential endpoints.")
  }

  # Init rates and check validity
  rate_exp <- x$para_exp[1]
  rate_ref <- x$para_ref[1]
  rate_pla <- x$para_pla[1]
  rates <- c(rate_exp, rate_ref, rate_pla)
  if (any(rates <= 0)) {
    stop("Rates must be positive.")
  }

  # Init probability for uncensored and check validity
  p_uncensor_exp <- x$para_exp[2]
  p_uncensor_ref <- x$para_ref[2]
  p_uncensor_pla <- x$para_pla[2]
  p_uncensor <- c(p_uncensor_exp, p_uncensor_ref, p_uncensor_pla)
  if (any((p_uncensor <= 0) | (p_uncensor > 1))) {
    stop("Probability for uncensored must be between zero and one.")
  }

  # Calculate effect
  effect <- sum(h(rates) * c(1, -Delta, Delta -1))
  if( effect > -1e-16 ){
    stop('Parameter vector is not located in the alternative.')
  }

  w <- allocation
  # Calculate variance
  var_teststat <- sum(c(1, Delta^2, (Delta -1)^2) / (w * p_uncensor))


  # Define 'method' for output
  method <- 'Power calculation for Wald-type test in three-arm trial with censored exponential endpoints'

  # Calculate missing parameter
  if (is.null(n)) {
    n <- ceiling( (qnorm(1-sig_level) + qnorm(power))^2 * var_teststat / effect^2)
    nExp <- round(n * w[1])
    nRef <- round(n * w[2])
    nPla <- round(n * w[3])
    n <- nExp + nRef + nPla
    if( any(c(nExp, nRef, nPla) / n != w) ){
      w <- c(nExp, nRef, nPla) / n
      note <- "'allocation' has been recalculated."
    }
    var_teststat <- sum(c(1, Delta^2, (Delta -1)^2) / (w * p_uncensor))
    power <- pnorm(qnorm(sig_level) - sqrt(n) * effect / sqrt(var_teststat), mean = 0, sd = 1)
  }
  if (is.null(power)) {
           power <- pnorm(qnorm(sig_level) - sqrt(n) * effect / sqrt(var_teststat), mean = 0, sd = 1)
  }
  if (is.null(sig_level)) {
          sig_level <- pnorm(qnorm(power) + sqrt(n) * effect / sqrt(var_teststat), mean = 0, sd = 1)
  }

  note <- NULL
  structure(list(method = method,
                 "Rate - Experiment" = rate_exp,
                 "Rate - Reference" = rate_ref,
                 "Rate - Placebo" = rate_pla,
                 "Prob Uncensor - Exp" = p_uncensor_exp,
                 "Prob Uncensor - Ref" = p_uncensor_ref,
                 "Prob Uncensor - Pla" = p_uncensor_pla,
                 n = n,
                 sig.level = sig_level,
                 power = power,
                 Delta = Delta,
                 allocation = w,
                 nExp = n*w[1],
                 nRef = n*w[2],
                 nPla = n*w[3],
                 note = note),
            class = "power.htest")

}
