#' @title estimate
#' @description generic function for calcuation the variance
#' @param x list with entries experiment, reference, placebo
#' @keywords internal
estimate <- function(x, Delta, ...) UseMethod("estimate")

#' @title calc_test_ret
#' @description generic function for performing the Wald-type test
#' @param x object of class negbin, poisson, exponential, ...
#' @keywords internal
calc_test_ret <- function(x, Delta, data_name, ...) UseMethod("calc_test_ret")

#' @title calc_power_ret
#' @description generic function for power related calculations for
#' the Wald-type test in three-arm trials
#' @param x object of class negbin, poisson, exponential, binary...
#' @keywords internal
calc_power_ret <- function(x, Delta, allocation, n, power, sig_level, distribution, ...) UseMethod("calc_power_ret")
