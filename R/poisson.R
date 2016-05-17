#' @title estimate.poisson
#' @description Parameter estimation for Wald-type test with Poisson distributed endpoints
#' @keywords internal
estimate.poisson <- function(x, Delta, ...) {

  var_estimation <- list(...)$var_estimation

  para <- numeric(3)
  para[1] <- mean(x$data_exp)
  para[2] <- mean(x$data_ref)
  para[3] <- mean(x$data_pla)

  effect <- sum(para * c(1, -Delta, Delta - 1))

  # if effect is in the null hypothesis, RML=ML
  # ML is faster and easier to calculate
  if (!is.null(var_estimation) && (var_estimation == "RML") && (effect >= 0)) {
    var_estimation <- "ML"
  }

  if (!is.null(var_estimation) && (var_estimation == "RML")) {

    # initilizing group sums and sample sizes
    n_exp <- x$n_exp
    n_ref <- x$n_ref
    n_pla <- x$n_pla
    sum_exp <- sum(x$data_exp)
    sum_ref <- sum(x$data_ref)
    sum_pla <- sum(x$data_pla)

    # Initialize matrix which defined restrictions to null hypothesis
    # and rates larger at least zero
    ui <- matrix(data = c(1, -Delta, (Delta-1),
                          diag(1, nrow = 3, ncol = 3)),
                 byrow = T,
                 ncol = 3)

    # starting value
    theta <- c(sum_exp/n_exp,
               Delta*sum_exp/n_exp,
               sum_exp/n_exp*(1-Delta^2)*Delta/(1-Delta))

    # Estimate rates restricted to null hypothesis
    COptimPar <- constrOptim(theta = theta,
                             f = poisson_likelihood,
                             grad = grad_poisson_likelihood,
                             ui = ui,
                             ci = c(0, 0, 0, 0),
                             n_exp = n_exp,
                             n_ref = n_ref,
                             n_pla = n_pla,
                             sum_exp = sum_exp,
                             sum_ref = sum_ref,
                             sum_pla = sum_pla)$par
    # output
    var_exp <- COptimPar[1]
    var_ref <- COptimPar[2]
    var_pla <- COptimPar[3]
    teststat_var <- x$n * (var_exp / x$n_exp +
                             Delta^2 * var_ref / x$n_ref +
                             (1-Delta)^2 * var_pla / x$n_pla)

  } else {
    teststat_var <- x$n * (para[1] / x$n_exp +
                             Delta^2 * para[2] / x$n_ref +
                             (1-Delta)^2 * para[3] / x$n_pla)
  }

  list(effect = effect,
       group_response = para,
       teststat_var = teststat_var)
}


#' @title calc_power_ret.poisson
#' @description Power related calculations for Wald-type test with Poisson distributed endpoints
#' @keywords internal
calc_power_ret.poisson <- function(x, Delta, allocation, n = NULL, power = NULL, sig_level = NULL, ...) {

  arguments <- list(...)

  var_estimation <- arguments$var_estimation
  if (is.null(var_estimation)) {
    var_estimation <- "ML"
  }

  if (any(c(length(x$para_exp), length(x$para_ref), length(x$para_pla)) != 1)) {
    stop("Only one parameter must be defined for power related calculations for poisson endpoints.")
  }

  rate_exp <- x$para_exp[1]
  rate_ref <- x$para_ref[1]
  rate_pla <- x$para_pla[1]
  rates <- c(rate_exp, rate_ref, rate_pla)
  if (any(rates <= 0)) {
    stop("Rates must be positive.")
  }

  # Calculate effect
  effect <- sum(rates * c(1, -Delta, Delta -1))
  if( effect > -1e-16 ){
    stop('Parameter vector is not located in the alternative.')
  }

  w <- allocation
  # Calculate variance
  var_teststat <- sum(rates * c(1, Delta^2, (Delta -1)^2) / w)
  var_teststat_rml <- limit_RML_poisson(rateExp1 = rate_exp,
                                        rateRef1 = rate_ref,
                                        ratePla1 = rate_pla,
                                        Delta = Delta,
                                        allocation = w)$sigma2_rml
  # Define 'method' for output
  method <- 'Power calculation for Wald-type test in three-arm trial with poisson endpoints'

  # Define 'method' for output
  switch(var_estimation,
         ML = ( method <- 'Power calculation for Wald-type test (with unrestriced variance estimation) in three-arm trial with Poisson endpoints' ),
         RML = ( method <- 'Power calculation for Wald-type test (with restriced variance estimation) in three-arm trial with Poisson endpoints' )
  )

  # Calculate missing parameter
  if (is.null(n)) {
    switch(var_estimation,
           ML =  (n <- ceiling( (qnorm(1-sig_level) + qnorm(power))^2 * var_teststat / effect^2) ),
           RML = (n <- ceiling( (qnorm(1-sig_level)*sqrt(var_teststat_rml)/ sqrt(var_teststat) + qnorm(power))^2 * var_teststat / effect^2) )
    )
    nExp <- round(n * w[1])
    nRef <- round(n * w[2])
    nPla <- round(n * w[3])
    n <- nExp + nRef + nPla
    if( any(c(nExp, nRef, nPla) / n != w) ){
      w <- c(nExp, nRef, nPla) / n
      note <- "'allocation' has been recalculated."
    }
    var_teststat <- sum(rates * c(1, Delta^2, (Delta -1)^2) / w)
    var_teststat_rml <- limit_RML_poisson(rateExp1 = rate_exp,
                                          rateRef1 = rate_ref,
                                          ratePla1 = rate_pla,
                                          Delta = Delta,
                                          allocation = w)$sigma2_rml
    switch(var_estimation,
           ML =  (power <- pnorm(qnorm(sig_level) - sqrt(n) * effect / sqrt(var_teststat), mean = 0, sd = 1) ),
           RML = (power <- pnorm((qnorm(sig_level)*sqrt(var_teststat_rml) - sqrt(n) * effect) / sqrt(var_teststat), mean = 0, sd = 1))
    )
  }
  if (is.null(power)) {
    switch(var_estimation,
           ML =  (power <- pnorm(qnorm(sig_level) - sqrt(n) * effect / sqrt(var_teststat), mean = 0, sd = 1) ),
           RML = (power <- pnorm((qnorm(sig_level)*sqrt(var_teststat_rml) - sqrt(n) * effect) / sqrt(var_teststat), mean = 0, sd = 1))
    )
  }
  if (is.null(sig_level)) {
    switch(var_estimation,
           ML =  (sig_level <- pnorm(qnorm(power) + sqrt(n) * effect / sqrt(var_teststat), mean = 0, sd = 1) ),
           RML = (sig_level <- pnorm((qnorm(power)*sqrt(var_teststat) + sqrt(n) * effect) / sqrt(var_teststat_rml), mean = 0, sd = 1))
    )
  }

  note <- NULL
  structure(list(method = method,
                 "Rate - Experiment" = rate_exp,
                 "Rate - Reference" = rate_ref,
                 "Rate - Placebo" = rate_pla,
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
