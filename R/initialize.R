#' @title init_test
#' @description Initialize object with the data from the three-arm trial which
#' will be used during the evaluation
#' @keywords internal
init_test <- function(data_exp,
                      data_ref,
                      data_pla,
                      distribution = NULL) {



  if (is.vector(data_exp) && is.vector(data_ref) && is.vector(data_pla)) {
    data_exp <- data_exp[!is.na(data_exp)]
    data_ref <- data_ref[!is.na(data_ref)]
    data_pla <- data_pla[!is.na(data_pla)]
    if (!is.numeric(c(data_exp, data_ref, data_pla))) {
      stop("Data must be numeric.")
    }
  } else if (is.matrix(data_exp) && is.matrix(data_ref) && is.matrix(data_pla)) {
    data_exp <- data_exp[!is.na(data_exp[, 1]) & !is.na(data_exp[, 2]),]
    data_ref <- data_ref[!is.na(data_ref[, 1]) & !is.na(data_ref[, 2]),]
    data_pla <- data_pla[!is.na(data_pla[, 1]) & !is.na(data_pla[, 2]),]
  } else {
    stop("Data not having the same type.")
  }

  x <- list(data_exp = data_exp,
            data_ref = data_ref,
            data_pla = data_pla)

  # Sample size and allocation
  n_exp <- length(data_exp)
  n_ref <- length(data_ref)
  n_pla <- length(data_pla)

  if (any(c(n_exp, n_ref, n_pla) < 2)) {
    stop("Each group must have at least two observations.", call. = FALSE)
  }

  x$n <- n_exp + n_ref + n_pla
  x$n_exp <- n_exp
  x$n_ref <- n_ref
  x$n_pla <- n_pla

  class(x) <- distribution
  x
}


#' @title init_power
#' @description Initialize object with the values for the power calculation
#' for a three-arm trial
#' @keywords internal
init_power <- function(para_exp,
                       para_ref,
                       para_pla,
                       distribution = NULL) {

  para_exp <- para_exp[!is.na(para_exp)]
  para_ref <- para_ref[!is.na(para_ref)]
  para_pla <- para_pla[!is.na(para_pla)]

  if (!is.numeric(c(para_exp, para_ref, para_pla))) {
    stop("The parameters must be numeric.")
  }
  length_inputs <- c(length(para_exp), length(para_ref), length(para_pla))

  if (any(length(para_exp) != length_inputs)) {
    stop("All distribution parameters must have the same length.")
  }


  x <- list(para_exp = para_exp,
            para_ref = para_ref,
            para_pla = para_pla)

  class(x) <- distribution
  x
}


