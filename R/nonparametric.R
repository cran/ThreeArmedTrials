#' @title estimate.nonparametric
#' @description Parameter estimation for Wald-type test in non-parametric model
#' @keywords internal
estimate.nonparametric <- function(x, Delta, ...) {

  estimate.normal(x = x, Delta = Delta, var_equal = FALSE)

}


#' @title calc_test_ret.nonparametric
#' @description Method for performing studentized permutation test of retention of effect hypothesis
#' in non-parametric model
#' @keywords internal
calc_test_ret.nonparametric <- function(x, Delta, data_name, ...) {

  n_perm <- list(...)$n_perm

  if (is.null(n_perm) || !is.naturalnumber(n_perm)) {
    stop('Number of permutations \'n_perm\' must be a natural number.')
  }

  if ((x$effect == 0) && (x$teststat_var == 0)) {
    warning("Effect and variance are zero. No p-value will be calculated.")
  }

  xExp <- x$data_exp
  xRef <- x$data_ref
  xPla <- x$data_pla

  # Sample size allocation
  nExp <- x$n_exp
  nRef <- x$n_ref
  nPla <- x$n_pla
  n <- x$n

  # Test statistic
  teststat <- sqrt(x$n) * x$effect / sqrt(x$teststat_var)


  # Permutation test
  obs <- c(xExp, xRef, xPla)

  xPerm <- t(sapply(X = 1:n_perm, function(i){sample(obs, replace = F)} ))
  xExpPerm <- xPerm[, 1:nExp]
  xRefPerm <- xPerm[, (nExp+1):(nExp+nRef)]
  xPlaPerm <- xPerm[, (nExp+nRef+1):n]

  xExpPermMean <- rowMeans(xExpPerm)
  xRefPermMean <- rowMeans(xRefPerm)
  xPlaPermMean <- rowMeans(xPlaPerm)

  sigma2ExpEst <- ( rowSums(xExpPerm^2) - nExp * xExpPermMean^2 ) / (nExp - 1)
  sigma2RefEst <- ( rowSums(xRefPerm^2) - nRef * xRefPermMean^2 ) / (nRef - 1)
  sigma2PlaEst <- ( rowSums(xPlaPerm^2) - nPla * xPlaPermMean^2 ) / (nPla - 1)

  sigma2_Tperm <- n * (sigma2ExpEst / nExp +
    Delta^2 * sigma2RefEst / nExp +
    (1-Delta)^2 * sigma2PlaEst / nExp )
  teststat_perm <- sqrt(n) *
    ( xExpPermMean - Delta * xRefPermMean + (Delta-1) * xPlaPermMean ) /
    sqrt(sigma2_Tperm)




  # Test statistic and p-value
  pvalue <- sum( teststat >= teststat_perm ) / (n_perm + 1) + 1 / (n_perm + 1)
  names(teststat) <- c('T')
  # Rename parameter vector
  names(x$group_response) <- c('Mean Exp', 'Mean Ref', 'Mean Pla')


  method_text <- "Studentized permutation test for the retention of effect hypothesis"


  structure(list(statistic = teststat,
                 p.value = pvalue,
                 method = method_text,
                 data.name = data_name,
                 estimate = x$group_response,
                 sample.size = x$n),
            class = "htest")

}


