#' @title estimate.negbin
#' @description Variance calculation for Wald-type test with negative binomial distributed endpoints
#' @keywords internal
estimate.negbin <- function(x, Delta, ...) {

  var_estimation <- list(...)$var_estimation

  para <- numeric(3)

  rate_exp <- para[1] <- mean(x$data_exp)
  rate_ref <- para[2] <- mean(x$data_ref)
  rate_pla <- para[3] <- mean(x$data_pla)

  effect <- sum(para * c(1, -Delta, Delta - 1))

  # if effect is in the null hypothesis, RML=ML
  # ML is faster and easier to calculate
  if (!is.null(var_estimation) && (var_estimation == "RML") && (effect >= 0)) {
    var_estimation <- "ML"
  }

  if (!is.null(var_estimation) && (var_estimation == "RML")) {

    # Initialize matrix which defined restrictions to null hypothesis
    # and rates larger at least zero
    ui <- matrix( c(1, -Delta, (Delta-1), 0,
                    diag(1, nrow=4, ncol = 4)),
                  byrow = T,
                  ncol = 4)

    # Starting value which is within the interior of the null hypothesis
    theta <- c(mean(x$data_exp),
               Delta * mean(x$data_exp),
               mean(x$data_exp) * (1-Delta^2) * Delta / (1-Delta),
               1)

    # Estimate rates restricted to null hypothesis
    opt_out <- constrOptim(theta = theta,
                           f = negbin_likelihood,
                           grad = NULL,
                           ui = ui,
                           ci = c(0, 0, 0, 0, 0),
                           x_exp = x$data_exp,
                           x_ref = x$data_ref,
                           x_pla = x$data_pla,
                           control = list(reltol = 1e-13,
                                          fnscale = -1))

    opt_count <- 1
    while ((opt_count <= 20) && (opt_out$convergence == 1)) {
      theta <- c(mean(x$data_exp),
                 Delta * mean(x$data_exp),
                 mean(x$data_exp) * (1-Delta^2) * Delta / (1-Delta),
                 runif(n = 1, min = 0.001, max = 20))

      # Estimate rates restricted to null hypothesis
      opt_out <- constrOptim(theta = theta,
                             f = negbin_likelihood,
                             grad = NULL,
                             ui = ui,
                             ci = c(0, 0, 0, 0, 0),
                             x_exp = x$data_exp,
                             x_ref = x$data_ref,
                             x_pla = x$data_pla,
                             control = list(reltol = 1e-13,
                                            fnscale = -1))
      opt_count <- opt_count + 1
      if (opt_count == 21) stop("RML cannot be calculated.")
    }

    var_exp <- opt_out$par[1] * (1 + opt_out$par[1] * opt_out$par[4])
    var_ref <- opt_out$par[2] * (1 + opt_out$par[2] * opt_out$par[4])
    var_pla <- opt_out$par[3] * (1 + opt_out$par[3] * opt_out$par[4])

  } else {

    n_exp <- x$n_exp
    n_ref <- x$n_ref
    n_pla <- x$n_pla
    obs <- c(x$data_exp, x$data_ref, x$data_pla)

    shape_ml <- 1/ .C("newton_Shape",
                      zufallszahlen=as.integer(obs),
                      mean1 = as.double(rate_exp),
                      mean2 = as.double(rate_ref),
                      mean3 = as.double(rate_pla),
                      n1 = as.integer(n_exp),
                      n2 = as.integer(n_ref),
                      n3 = as.integer(n_pla),
                      theta_out = as.double(numeric(1)) )$theta_out
    var_exp <- rate_exp * (1 + rate_exp * shape_ml)
    var_ref <- rate_ref * (1 + rate_ref * shape_ml)
    var_pla <- rate_pla * (1 + rate_pla * shape_ml)
  }

  teststat_var <- x$n * (var_exp / x$n_exp +
                           Delta^2 * var_ref / x$n_ref +
                           (1-Delta)^2 * var_pla / x$n_pla)

  list(effect = effect,
       group_response = para,
       teststat_var = teststat_var)
}




#' @title calc_power_ret.negbin
#' @description Power related calculations for Wald-type test with Poisson distributed endpoints
#' @keywords internal
calc_power_ret.negbin <- function(x, Delta, allocation, n = NULL, power = NULL, sig_level = NULL, ...) {

  arguments <- list(...)

  var_estimation <- arguments$var_estimation
  if (is.null(var_estimation)) {
    var_estimation <- "ML"
  }

  if (any(c(length(x$para_exp), length(x$para_ref), length(x$para_pla)) != 2)) {
    stop("Two parameters must be defined for power related calculations for negative binomial endpoints.")
  }

  rate_exp <- x$para_exp[1]
  rate_ref <- x$para_ref[1]
  rate_pla <- x$para_pla[1]
  rates <- c(rate_exp, rate_ref, rate_pla)
  if (any(rates <= 0)) {
    stop("Rates must be positive.")
  }

  shape_exp <- x$para_exp[2]
  shape_ref <- x$para_ref[2]
  shape_pla <- x$para_pla[2]
  shapes <- c(shape_exp, shape_ref, shape_pla)
  if (any((shapes <= 0) || (shape_exp != shapes))) {
    stop("Shape parameters must be positive and identical.")
  }
  shape <- unique(shapes)

  # Calculate effect
  effect <- sum(rates * c(1, -Delta, Delta -1))
  if( effect > -1e-16 ){
    stop('Parameter vector is not located in the alternative.')
  }

  w <- allocation
  # Calculate variance
  var_teststat <- sum( rates * (1 + rates * shape) * c(1, Delta^2, (1-Delta)^2) / w )
  if (var_estimation == "RML") {
  var_teststat_rml <- taNegbin.LimitRestMLE(rateExp1 = rate_exp,
                                            rateRef1 = rate_ref,
                                            ratePla1 = rate_pla,
                                            shape1 = shape,
                                            Delta = Delta,
                                            allocation = w)$sigma2.rest
  }

  # Define 'method' for output
  switch(var_estimation,
         ML = ( method <- 'Power calculation for Wald-type test (with unrestriced variance estimation) in three-arm trial with negative binomial endpoints' ),
         RML = ( method <- 'Power calculation for Wald-type test (with restriced variance estimation) in three-arm trial with negative binomial endpoints' )
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
    var_teststat <- sum( rates * (1 + rates * shape) * c(1, Delta^2, (1-Delta)^2) / w )
    var_teststat_rml <- taNegbin.LimitRestMLE(rateExp1 = rate_exp,
                                              rateRef1 = rate_ref,
                                              ratePla1 = rate_pla,
                                              shape1 = shape,
                                              Delta = Delta,
                                              allocation = w)$sigma2.rest
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
                 shape = shape,
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

