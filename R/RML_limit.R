#' @title limit_RML_binary
#' @description Calculate limit of RML for Poisson endpoints
#' @param pExp1 Experimental treatment group rate in the (assumed) alternative
#' @param pRef1 Experimental treatment group rate in the (assumed) alternative
#' @param pPla1 Experimental treatment group rate in the (assumed) alternative
#' @param Delta Non-inferiority/Superiority margin
#' @param allocation Sample size allocation
#' @param h Function used when defining hypothesis
#' @param h_inv Inverse function to h
#' @keywords internal
limit_RML_binary <- function(pExp1, pRef1, pPla1, Delta, allocation = c(1/3, 1/3, 1/3), h, h_inv){

  if (is.null(h) || is.null(h_inv)) {
    h <- identity
    h_inv <- identity
  }

  # M. Mielke "Maximum Likelihood Theory for Retention of Effect Non-inferiority Trials"
  ## Page 34
  kl_divergence <- function(paraH0, pH1, Delta, w){
    pH0 <- c(h_inv(Delta * h(paraH0[1]) + (1 - Delta) * h(paraH0[2])), paraH0[1], paraH0[2])
    if (any((pH0 <=0) | (pH0 >= 1))) {
      return(99999)
    }
    sum(w * (pH1 * log(pH1 / pH0) + (1-pH1) * log((1-pH1) / (1-pH0))))
  }
  opt_out <- optim(par = c(pRef1, pPla1),
                   kl_divergence,
                   pH1 = c(pExp1, pRef1, pPla1),
                   w = allocation,
                   Delta = Delta,
                   method = "L-BFGS-B",
                   lower = c(0, 0),
                   upper = c(1,1))
  conv_counter <- 1
  while ((opt_out$convergence != 0) && (conv_counter <= 5)) {
    opt_out <- optim(par = runif(2),
                     kl_divergence,
                     pH1 = c(pExp1, pRef1, pPla1),
                     w = allocation,
                     Delta = Delta,
                     method = "L-BFGS-B",
                     lower = c(0, 0),
                     upper = c(1,1))
    conv_counter <- conv_counter + 1
  }
  if (conv_counter > 5) {
    stop("Minimization of Kullback-Leibler did not convergence.")
  }

  pRef0 <- opt_out$par[1]
  pPla0 <- opt_out$par[2]
  pExp0 <- h_inv(Delta * h(pRef0) + (1 - Delta) * h(pPla0))
  p <- c(pExp0, pRef0, pPla0)
  group_var <- numDeriv::grad(func = h, x = p)^2 * p * (1-p)
  sigma2_rml <- sum(group_var * c(1, Delta^2, (1-Delta)^2) / allocation)

  list(pExp0 = pExp0,
       pRef0 = pRef0,
       pPla0 = pPla0,
       sigma2_rml = sigma2_rml)
}



#' @title limit_RML_poisson
#' @description Calculate limit of RML for Poisson endpoints
#' @param rateExp1 Experimental treatment group rate in the (assumed) alternative
#' @param rateRef1 Experimental treatment group rate in the (assumed) alternative
#' @param ratePla1 Experimental treatment group rate in the (assumed) alternative
#' @param Delta Non-inferiority/Superiority margin
#' @param allocation Sample size allocation
#' @keywords internal
limit_RML_poisson <- function(rateExp1, rateRef1, ratePla1, Delta, allocation = c(1/3, 1/3, 1/3)){


  kl_divergence <- function(paraH0, ratesH1, Delta, w){
    ratesH0 <- c(Delta * paraH0[1] + (1 - Delta) * paraH0[2], paraH0[1], paraH0[2])
    sum(w * (ratesH0 - ratesH1 + ratesH1 * (log(ratesH1) - log(ratesH0))))
  }

  opt_out <- optim(par = c(rateRef1, ratePla1),
                   kl_divergence,
                   ratesH1 = c(rateExp1, rateRef1, ratePla1),
                   w = allocation,
                   Delta = Delta,
                   method = "L-BFGS-B",
                   lower = c(0, 0))
  if (opt_out$convergence != 0) {
    stop("Minimization of Kullback-Leibler did not convergence.")
  }

  rateRef0 <- opt_out$par[1]
  ratePla0 <- opt_out$par[2]
  rateExp0 <- Delta * rateRef0 + (1 - Delta) * ratePla0

  sigma2_rml <- rateExp0 / allocation[1] +
    Delta^2 * rateRef0 / allocation[2] +
    (1-Delta)^2 * ratePla0 / allocation[3]

  list(rateExp0 = rateExp0,
       rateRef0 = rateRef0,
       ratePla0 = ratePla0,
       sigma2_rml = sigma2_rml)
}



# title Limit of restricted maximum-likelihood estimator in case of negative binomial distributed endpoints
# param rateExp1 A numeric value specifying the rate of the experimental treatment group in the alternative hypothesis
# param rateRef1 A numeric value specifying the rate of the reference treatment group in the alternative hypothesis
# param ratePla1 A numeric value specifying the rate of the placebo treatment group in the alternative hypothesis
# param shape1 A numeric value specifying the shape parameter
# param Delta A numeric value specifying the non-inferiority/superiority margin
# param allocation A (non-empty) vector specifying the sample size allocation (wExp, wRef, wPla)
# return A list containing the following components:
# item{rateExp0, rateRef0, rateRef0}{The limit of the maximum-likelihood estimator for the rates when estimed restricted to the boundary of the null hypothesis}
# item{shape0}{The limit of the maximum-likelihood estimator for the shape parameter when estimed restricted to the boundary of the null hypothesis}
# item{sigma2.rest}{The limit of the maximum-likelihood variance estimator for the Wald-type test when restricted to the boundary of the null hypothesis}
taNegbin.LimitRestMLE <- function(rateExp1, rateRef1, ratePla1, shape1, Delta, allocation = c(1/3, 1/3, 1/3)){

  KL.Divergenz <- function(zeta, para.alternat, Delta, w){

    rateExp1 <- para.alternat[1]
    rateRef1 <- para.alternat[2]
    ratePla1 <- para.alternat[3]
    shape1 <- para.alternat[4]

    xLimit <- 1
    while(pnbinom(q=xLimit, mu=max(ratePla1, rateExp1, rateRef1), size = 1/max(shape1,zeta[3])) < 1){
      xLimit <- xLimit * 5
    }

    x <- seq(from=0, to=xLimit)

    rateExp <- zeta[1]
    rateRef <- zeta[2]
    shape <- zeta[3]
    if( rateExp <= 0 || rateRef <= 0 || shape <= 0){
      return(Inf)
    }
    ratePla <- ( rateExp - Delta * rateRef ) / (1-Delta)
    if(ratePla <= 0){
      return(Inf)
    }
    KL.Exp <- sum( ( dnbinom(x, mu = rateExp1, size = 1 / shape1, log = TRUE) -
                       dnbinom(x, mu = rateExp, size = 1 / shape, log = TRUE) ) *
                     dnbinom(x, mu = rateExp1, size = 1 / shape1, log = FALSE) )
    KL.Ref <- sum( ( dnbinom(x, mu = rateRef1, size = 1 / shape1, log = TRUE) -
                       dnbinom(x, mu = rateRef, size = 1 / shape, log = TRUE) ) *
                     dnbinom(x, mu = rateRef1, size = 1 / shape1, log = FALSE) )
    KL.Pla <- sum( ( dnbinom(x, mu = ratePla1, size = 1 / shape1, log = TRUE) -
                       dnbinom(x, mu = ratePla, size = 1 / shape, log = TRUE) ) *
                     dnbinom(x, mu = ratePla1, size = 1 / shape1, log = FALSE) )

    return( w[1]*KL.Exp + w[2]*KL.Ref + w[3]*KL.Pla )
  }

  opt.KL <- optim(c(rateRef1, rateRef1, shape1), KL.Divergenz,
                  para.alternat = c(rateExp1, rateRef1, ratePla1, shape1),
                  Delta = Delta,
                  w = allocation)$par

  rateExp0 <- opt.KL[1]
  rateRef0 <- opt.KL[2]
  ratePla0 <- ( rateExp0 - Delta * rateRef0 ) / ( 1 - Delta )
  shape0 <- opt.KL[3]
  sigma2.rest <- rateExp0 * ( 1 + rateExp0 * shape0 ) / allocation[1] + Delta^2 * rateRef0 * ( 1 + rateRef0 * shape0 ) / allocation[2] + (1-Delta)^2 * ratePla0 * ( 1 + ratePla0 * shape0 ) / allocation[3]
  return(list(
    rateExp0 = rateExp0,
    rateRef0 = rateRef0,
    ratePla0 = ratePla0,
    shape0 = shape0,
    sigma2.rest = sigma2.rest
  )
  )
}
