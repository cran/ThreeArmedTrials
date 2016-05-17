#' @title estimate.binary
#' @description Parameter estimation for Wald-type test with binary endpoints
#' @keywords internal
estimate.binary <- function(x, Delta, ...) {

  var_estimation <- list(...)$var_estimation
  h <- list(...)$h
  h_inv <- list(...)$h_inv

  # Check data
  if (any( !(c(x$data_exp, x$data_ref, x$data_pla) %in% c(0,1)))) {
    stop("Binary data must be either 0 or 1.")
  }

  # Set h to identity if not defined
  if (is.null(h) || is.null(h_inv)) {
    h <- identity
    h_inv <- identity
    warning("Either h or h_inv is NULL.")
  }

  # initialize group specific sample sizes
  n_exp <- x$n_exp
  n_ref <- x$n_ref
  n_pla <- x$n_pla

  para <- numeric(3)
  para[1] <- mean(x$data_exp)
  para[2] <- mean(x$data_ref)
  para[3] <- mean(x$data_pla)

  effect <- h(para[1]) - Delta * h(para[2]) + (Delta - 1) * h(para[3])

  # if effect is in the null hypothesis, RML=ML
  # ML is faster and easier to calculate
  if (!is.null(var_estimation) && (var_estimation == "RML") && (effect >= 0)) {
    var_estimation <- "ML"
  }

  if (!is.null(var_estimation) && (var_estimation == "RML")) {

    # initilizing group sums and sample sizes
    sum_exp <- sum(x$data_exp)
    sum_ref <- sum(x$data_ref)
    sum_pla <- sum(x$data_pla)

    COptimPar <- optim(par = c(mean(x$data_ref), mean(x$data_pla)),
                       fn = binary_likelihood_rest,
                       Delta = Delta,
                       n_exp = n_exp,
                       n_ref = n_ref,
                       n_pla = n_pla,
                       sum_exp = sum_exp,
                       sum_ref = sum_ref,
                       sum_pla = sum_pla,
                       h = h,
                       h_inv = h_inv,
                       method = "L-BFGS-B",
                       lower = c(0,0),
                       upper = c(1,1))$par
    rest_para <- c(h_inv(Delta * h(COptimPar[1]) + (1-Delta) * h(COptimPar[2])), COptimPar[1:2])
    group_var <- numDeriv::grad(func = h, x = rest_para)^2 * rest_para * (1-rest_para)
  } else {
    group_var <- numDeriv::grad(func = h, x = para)^2 * para * (1-para)
  }

  teststat_var <- x$n * sum(group_var * c(1, Delta^2, (1-Delta)^2) * c(1/n_exp, 1/n_ref, 1/n_pla))

  list(effect = effect,
       group_response = para,
       teststat_var = teststat_var)
}


#' @title calc_power_ret.binary
#' @description Power related calculations for Wald-type test with binary endpoints
#' @keywords internal
calc_power_ret.binary <- function(x, Delta, allocation, n = NULL, power = NULL, sig_level = NULL, ...) {

  arguments <- list(...)
  var_estimation <- arguments$var_estimation
  h <- list(...)$h
  h_inv <- list(...)$h_inv

  if (is.null(var_estimation)) {
    var_estimation <- "ML"
  }
  # Set h to identity if not defined
  if (is.null(h) || is.null(h_inv)) {
    h <- identity
    h_inv <- identity
    warning("Either h or h_inv is NULL.")
  }

  if (any(c(length(x$para_exp), length(x$para_ref), length(x$para_pla)) != 1)) {
    stop("Only one parameter must be defined for power related calculations for binary endpoints.")
  }

  p_exp <- x$para_exp[1]
  p_ref <- x$para_ref[1]
  p_pla <- x$para_pla[1]
  p <- c(p_exp, p_ref, p_pla)
  if (any(p <= 0)) {
    stop("Probabilities must be positive.")
  }

  # Calculate effect
  effect <- h(p_exp) - Delta * h(p_ref) + (Delta - 1) * h(p_pla)
  if( effect > -1e-16 ){
    stop('Parameter vector is not located in the alternative.')
  }

  w <- allocation
  # Calculate variance
  group_var <- numDeriv::grad(func = h, x = p)^2 * p * (1-p)
  var_teststat <- sum(group_var * c(1, Delta^2, (1-Delta)^2) / w)
  var_teststat_rml <- limit_RML_binary(pExp1 = p_exp,
                                       pRef1 = p_ref,
                                       pPla1 = p_pla,
                                       Delta = Delta,
                                       allocation = w,
                                       h = h,
                                       h_inv = h_inv)$sigma2_rml
  # Define 'method' for output
  switch(var_estimation,
         ML = ( method <- 'Power calculation for Wald-type test (with unrestriced variance estimation) in three-arm trial with binary endpoints' ),
         RML = ( method <- 'Power calculation for Wald-type test (with restriced variance estimation) in three-arm trial with binary endpoints' )
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
    var_teststat <- sum(group_var * c(1, Delta^2, (1-Delta)^2) / w)
    var_teststat_rml <- limit_RML_binary(pExp1 = p_exp,
                                         pRef1 = p_ref,
                                         pPla1 = p_pla,
                                         Delta = Delta,
                                         allocation = w,
                                         h = h,
                                         h_inv = h_inv)$sigma2_rml

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
                 "Rate - Experimental" = p_exp,
                 "Rate - Reference" = p_ref,
                 "Rate - Placebo" = p_pla,
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
