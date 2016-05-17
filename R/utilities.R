#' @title check_missing
#' @description Check if all arguments are defined
#' @param args Character vector of arguments to be checked for existence.
#' @param envir Environment in which the arguments are defined.
check_missing <- function(args = NULL, envir = parent.frame()) {

  # Check which are missing
  isMissing <- sapply(args, function(arg){ eval(call("missing", as.name(arg)), envir = envir) })

  # Return error message when missing is not fulfilled
  if (any(isMissing)) {
    args_missing <- args[isMissing]

    if (length(args_missing) == 1) {
      stop( paste("Argument", args_missing, "is missing.", sep = " "),
            call. = FALSE)
    } else if (length(args_missing) > 1) {
      stop( paste("Arguments", paste(args_missing, collapse = ', '),
                  "are missing.", sep = " "), call. = FALSE)
    } else{
      return(invisible())
    }
  }
}

#' @title check_RET_arguments
#' @description Check arguments for their respective condition
#' @param Delta A numeric value specifying the non-inferiority or superiority margin.
#' Is between 0 and 1 in case of non-inferiority and larger than 1 in case of superiority.
#' @param sig.level A numeric value specifying the significance level (type I error probability)
#' @param power A numeric value specifying the target power (1 - type II error probability)
#' @param n The total sample size. Needs to be at least 7.
#' @param allocation A (non-empty) vector specifying the sample size allocation (nExp/n, nRef/n, nPla/n)
check_RET_arguments <- function(sig.level, power, Delta, n, allocation) {

  if( !is.null(sig.level) && !is.numeric(sig.level) || any(0 >= sig.level | sig.level >= 1) ){
    stop("Significance level 'sig.level' must be numeric in (0, 1).")
  }

  if( !is.null(power) && !is.numeric(power) || any(0 >= power | power >= 1) ){
    stop("'power' must be numeric in (0, 1).")
  }

  if( !is.null(Delta) && !is.numeric(Delta) || (0 >= Delta) ){
    stop("'Delta' must be larger than 0.")
  }

  if( !is.null(n) && !is.numeric(n) || any(6 >= n) ){
    stop("'n' must be larger than 6.")
  }

  w <- allocation
  if (!is.null(w) && !is.numeric(w) || any(w >= 1, w <= 0, abs(sum(w)-1) > 1e-10, length(w) != 3)) {
    stop("'allocation' must not have length 3, sum up to 1, and have only entries between 0 and 1.")
  }
}

#' @title loglikelihood_binary
#' @description log likelihood of Bernoulli function
#' @param p numeric vector of probabilities with length 3
#' @param xExp numeric vector of probabilities with length 3
#' @param xRef numeric vector of probabilities with length 3
#' @param xPla numeric vector of probabilities with length 3
loglikelihood_binary <- function(p, xExp, xRef, xPla) {
  pExp <- p[1]
  pRef <- p[2]
  pPla <- p[3]

  nExp <- length(xExp)
  nRef <- length(xRef)
  nPla <- length(xPla)

  sumExp <- sum(xExp)
  sumRef <- sum(xRef)
  sumPla <- sum(xPla)

  log_l <- sumExp * log(pExp) + (nExp - sumExp) * log(1 - pExp) +
    sumRef * log(pRef) + (nRef - sumRef) * log(1 - pRef) +
    sumPla * log(pPla) + (nPla - sumPla) * log(1 - pPla)

  -log_l
}

#' @title is.naturalnumber
#' @description check if input is natural number
#' @param x numeric number to be checked
#' @param tol maximum accepted tolerance when checking if natural
is.naturalnumber <- function(x, tol = .Machine$double.eps^0.5) {
  ((abs(x - round(x)) < tol) & (x > 0))
}
