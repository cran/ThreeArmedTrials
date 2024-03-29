% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/waldtype_optAlloc.R
\name{opt_alloc_RET}
\alias{opt_alloc_RET}
\title{Optimal sample size for three-arm trials when analyzed with a Wald-type test}
\usage{
opt_alloc_RET(experiment, reference, placebo, Delta, distribution, h = NULL)
}
\arguments{
\item{experiment}{a numeric vector specifying the parameters of the experimental
treatment group in the alternative hypothesis}

\item{reference}{a numeric vector specifying the parameters of the reference
treatment group in the alternative hypothesis}

\item{placebo}{a numeric vector specifying the parameters of the placebo
treatment group in the alternative hypothesis}

\item{Delta}{a numeric value specifying the non-inferiority/superiority margin}

\item{distribution}{a character specifying the distribution of the endpoints. Must
must be either of \code{"poisson"}, \code{"negbin"}, \code{"exponential"}, \code{"normal"}}

\item{h}{Function measuring the efficacy; used to defined hypothesis}
}
\value{
Vector with optimal sample size allocation in the order
(experiment, reference, placebo)
}
\description{
Calculate optimal sample size allocation for Wald-type test for
superiority or non-inferiority of the experimental treatment versus
reference treatment with respect to placebo
}
\details{
The arguments \code{experiment}, \code{reference}, and \code{placebo} define
the parameters of the endpoint distribution for the respective groups: \cr
\code{distribution = "poisson"}: \code{experiment}, \code{reference}, and
\code{placebo} must have length one and define the means. \cr
\code{distribution = "negbin"}: \code{experiment}, \code{reference}, and
\code{placebo} must have length two and define the mean in the first entry
and the shape parameter in the second entry. \cr
\code{distribution = "exponential"}: \code{experiment}, \code{reference}, and
\code{placebo} must have length two and define the mean in the first entry
and the probability for an uncensored observation in the second entry. \cr
\code{distribution = "normal"}: \code{experiment}, \code{reference}, and
\code{placebo} must have length two and define the mean in the first entry
and the variance in the second entry.
}
\examples{
opt_alloc_RET(experiment = 1,
              reference = 1,
              placebo = 3,
              Delta = 0.8,
              distribution = "poisson")

}
\keyword{allocation}
\keyword{waldtype}
