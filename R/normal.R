#' @title estimate.normal
#' @description Parameter estimation for Wald-type test with normally distributed endpoints
#' @keywords internal
estimate.normal <- function(x, Delta, ...) {

  var_equal <- list(...)$var_equal

  if (is.null(var_equal)) {
    warning("Parameter \'var_equal\' not specified. Heteroscedastic group variances are assumed.")
    var_equal <- FALSE
  }


  if (var_equal) {

    teststat_var <- x$n * (1 / x$n_exp + Delta^2 / x$n_ref + (1-Delta)^2 / x$n_pla) *
      ( (x$n_exp-1) * var(x$data_exp) +
          (x$n_ref-1) * var(x$data_ref) +
          (x$n_pla-1) * var(x$data_pla) ) /
      (x$n - 3)
    degF <- x$n - 3

    para <- numeric(3)
    para[1] <- mean(x$data_exp)
    para[2] <- mean(x$data_ref)
    para[3] <- mean(x$data_pla)

    effect <- sum(para * c(1, -Delta, Delta - 1))

    return(
      list(effect = effect,
           group_response = para,
           teststat_var = teststat_var)
    )

  } else {
    return(estimate.default(x = x, Delta = Delta))
  }

}


#' @title calc_test_ret.normal
#' @description Method for performing test of retention of effect hypothesis for
#' normally distributed endpoints
#' @keywords internal
calc_test_ret.normal <- function(x, Delta, data_name, ...) {

  var_equal <- list(...)$var_equal
  if (is.null(var_equal)) {
    var_equal <- FALSE
  }

  # Degrees of freedom
  # Variance estimator
  if (var_equal) {
    degF <- x$n - 3
  } else {
    if (x$teststat_var == 0) {
      degF <- 1
    } else {
      degF <- x$teststat_var^2 / x$n^2 /
        ( var(x$data_exp)^2 / (x$n_exp^2 * (x$n_exp - 1)) +
            Delta^4 * var(x$data_ref)^2 / (x$n_ref^2 * (x$n_ref - 1)) +
            (1-Delta)^4 * var(x$data_pla)^2 / (x$n_pla^2 * (x$n_pla - 1))  )
    }
  }

  if ((x$effect == 0) && (x$teststat_var == 0)) {
    warning("Effect and variance are zero. No p-value will be calculated.")
  }


  # Test statistic and p-value
  teststat <- sqrt(x$n) * x$effect / sqrt(x$teststat_var)
  pvalue <- pt(teststat, df = degF)
  names(teststat) <- c('T')
  names(degF) <- "df"
  # Rename parameter vector
  names(x$group_response) <- c('Mean Exp', 'Mean Ref', 'Mean Pla')

  if (var_equal) {
    method_text <- paste0("Wald-type test for the retention of effect hypothesis ",
                          "for normally distributed endpoints with equal variances.")
  } else {
    method_text <- paste0("Wald-type test for the retention of effect hypothesis ",
                          "for normally distributed endpoints with unequal variances.")
  }


  structure(list(statistic = teststat,
                 p.value = pvalue,
                 method = method_text,
                 data.name = data_name,
                 estimate = x$group_response,
                 parameter = degF,
                 sample.size = x$n),
            class = "htest")

}





calc_power_ret.normal <- function(x, Delta, allocation, n = NULL, power = NULL, sig_level = NULL, ...){

  arguments <- list(...)
  var_equal <- arguments$var_equal
  if (is.null(var_equal)) {
    var_equal <- FALSE
  }

  if (any(c(length(x$experiment), length(x$reference), length(x$placebo)) != 2)) {
    stop("Two parameters must be defined for power related calculations for normal endpoints.")
  }

  mean_exp <- x$para_exp[1]
  mean_ref <- x$para_ref[1]
  mean_pla <- x$para_pla[1]
  means <- c(mean_exp, mean_ref, mean_pla)

  var_exp <- x$para_exp[2]
  var_ref <- x$para_ref[2]
  var_pla <- x$para_pla[2]
  variances <- c(var_exp, var_ref, var_pla)

  if (any(variances <= 0)) {
    stop("Variances must be positive.")
  }

  # Calculate effect
  effect <- sum(means * c(1, -Delta, Delta -1))
  if( effect > -1e-16 ){
    stop('Parameter vector is not located in the alternative.')
  }

  w <- allocation
  # Calculate variance
  var_teststat <- sum(means * c(1, Delta^2, (Delta -1)^2) / w)


  tol <- .Machine$double.eps^0.25
  p.body <-  quote({
    n_exp <- round(n * w[1])
    n_ref <- round(n * w[2])
    n_pla <- round(n * w[3])
    ncp <- effect / sqrt(var_exp / n_exp +
                           Delta^2 * var_ref / n_ref +
                           (1-Delta)^2 * var_pla / n_pla)

    if (var_equal) {
      degF <- n_exp + n_ref + n_pla - 3
    } else {
      degF <- (var_exp / n_exp + Delta^2 * var_ref / n_ref + (1-Delta)^2 * var_pla / n_pla)^2 /
        ( var_exp^2 / (n_exp^2 * (n_exp-1)) +
            Delta^4 * var_ref^2 / (n_ref^2 * (n_ref-1)) +
            (1-Delta)^4 * var_pla^2 / (n_pla^2 * (n_pla-1))  )
    }

    critical_val <- qt(sig_level, df = degF, lower.tail = TRUE)
    pt(critical_val, df = degF, ncp = ncp, lower.tail = TRUE)
  })

  if (is.null(power)) {
    power <- eval(p.body)
  }  else if (is.null(n)) {
    n <- uniroot(function(n) eval(p.body) - power, c(ceiling(1/min(w))+3, 1e+07),
                 tol = tol, extendInt = "upX")$root
    n_exp <- round(n * w[1])
    n_ref <- round(n * w[2])
    n_pla <- round(n * w[3])
    n <- n_exp + n_ref + n_pla
    n_group <- c(n_exp, n_ref, n_pla)
    if( any(n_group / n != w) ){
      w <- n_group / n
      note <- "'allocation' has been recalculated."
    }
    power <- eval(p.body)
  } else if (is.null(sig_level)) {
    sig_level <- uniroot(function(sig_level) power - eval(p.body),
                         c(1e-4, 1 - 1e-4), tol = tol, extendInt = "yes")$root
  } else {
    stop("internal error", domain = NA)
  }

  # Define 'method' for output
  method <- 'Power calculation for Wald-type test in three-arm trial with poisson endpoints'
  note <- NULL
  structure(list(method = method,
                 meanExp = mean_exp,
                 meanRef = mean_ref,
                 meanPla = mean_pla,
                 varExp = var_exp,
                 varRef = var_ref,
                 varPla = var_pla,
                 n = n,
                 sig.level = sig_level,
                 power = power,
                 Delta = Delta,
                 allocation = w,
                 nExp = n * w[1],
                 nRef = n * w[2],
                 nPla = n * w[3],
                 note = note),
            class = "power.htest")
}

