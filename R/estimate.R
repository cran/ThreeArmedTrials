#' @title estimate
#' @description generic function for calculating the variance
#' @param x list with entries experiment, reference, placebo
#' @keywords internal
estimate <- function(x, Delta, ...) UseMethod("estimate")
