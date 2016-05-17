#' @title estimate.default
#' @description Default variance calculation for Wald-type test
#' @keywords internal
estimate.default <- function(x, Delta, ...) {

  para <- numeric(3)
  para[1] <- mean(x$data_exp)
  para[2] <- mean(x$data_ref)
  para[3] <- mean(x$data_pla)

  effect <- sum(para * c(1, -Delta, Delta - 1))

  var_exp <- var(x$data_exp)
  var_ref <- var(x$data_ref)
  var_pla <- var(x$data_pla)

  teststat_var <- x$n * (var_exp / x$n_exp +
                           Delta^2 * var_ref / x$n_ref +
                           (1-Delta)^2 * var_pla / x$n_pla)

  list(effect = effect,
       group_response = para,
       teststat_var = teststat_var)
}

#' @title test_ret.default
#' @description Default method for performing Wald-type test
#' @keywords internal
calc_test_ret.default <- function(x, Delta, data_name, ...) {

  var_estimation <- list(...)$var_estimation

  if ((x$effect == 0) && (x$teststat_var == 0)) {
    warning("Effect and variance are zero. No p-value will be calculated.")
  }

  # Test statistic and p-value
  teststat <- sqrt(x$n) * x$effect / sqrt(x$teststat_var)
  pvalue <- pnorm(teststat)
  names(teststat) <- c('T')

  # Rename parameter vector
  names(x$group_response) <- c('Mean Exp', 'Mean Ref', 'Mean Pla')

  if (class(x) != "list") {
    var_estimation <- list(...)$var_estimation
    if (!is.null(var_estimation)) {
      var_text <- switch(var_estimation,
                         ML = {' (with unrestriced variance estimation) '},
                         RML = {' (with restriced variance estimation) '})
    } else {
      var_text <- ""
    }
    method_text <- paste0("Wald-type test for the retention of effect hypothesis ",
                          var_text,
                          "for the \'",
                          class(x), "\ model'")
  } else {
    method_text <- 'Wald-type test for the retention of effect hypothesis using method of moments estimators'
  }

  structure(list(statistic = teststat,
                 p.value = pvalue,
                 method = method_text,
                 data.name = data_name,
                 estimate = x$group_response,
                 sample.size = x$n),
            class = "htest")
}


#' @title power_ret.default
#' @description Default method for power calculations in three-arm trials
#' @keywords internal
calc_power_ret.default <- function(x, Delta, allocation, n = NULL, power = NULL, sig_level = NULL, method = NULL, ...) {

  mean_exp <- x$para_exp[1]
  mean_ref <- x$para_ref[1]
  mean_pla <- x$para_pla[1]

  var_exp <- x$para_exp[2]
  var_ref <- x$para_ref[2]
  var_pla <- x$para_pla[2]

  if (any(c(var_exp, var_ref, var_pla) <= 0)) {
    stop("Variances must be positive.")
  }

  effect <- mean_exp - Delta * mean_ref + (Delta - 1) * mean_pla
  if( effect > -1e-16 ){
    stop('Parameter vector is not located in the alternative.')
  }

  # initialize note
  note <- NULL

  # initialize w
  w <- allocation

  # Calculate variance
  var_teststat <- var_exp / w[1] +
    Delta^2 * var_ref / w[2] +
    (1-Delta)^2 * var_pla / w[3]

  # Define 'method' for output
  if (is.null(method)) {
    method <- 'Wald-type test for three-arm trials power calculation'
  }

  # Calculate missing parameter
  if (is.null(n)) {
    n <- ceiling( (qnorm(1-sig_level) + qnorm(power))^2 * var_teststat / effect^2 )
    nExp <- round(n * w[1])
    nRef <- round(n * w[2])
    nPla <- round(n * w[3])
    n <- nExp + nRef + nPla
    if( any(c(nExp, nRef, nPla) / n != w) ){
      w <- c(nExp, nRef, nPla) / n
      note <- "'allocation' was adjusted to match the total sample size."
    }
    var_teststat <- var_exp / w[1] +
      Delta^2 * var_ref / w[2] +
      (1-Delta)^2 * var_pla / w[3]
    power <- pnorm(qnorm(sig_level) - sqrt(n) * effect / sqrt(var_teststat), mean = 0, sd = 1)
  }
  if (is.null(power)) {
    power <- pnorm(qnorm(sig_level) - sqrt(n) * effect / sqrt(var_teststat), mean = 0, sd = 1)
  }
  if (is.null(sig_level)) {
    sig.level <- pnorm(qnorm(power) + sqrt(n) * effect / sqrt(var_teststat), mean = 0, sd = 1)
  }

  structure(list(method = method,
                 "Expected Value - Experiment" = mean_exp,
                 "Expected Value - Reference" = mean_ref,
                 "Expected Value - Placebo" = mean_pla,
                 "Variance - Experiment" = var_exp,
                 "Variance - Reference" = var_ref,
                 "Variance - Placebo" = var_pla,
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
