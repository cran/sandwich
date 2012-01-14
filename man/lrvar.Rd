\name{lrvar}
\alias{lrvar}

\title{Long-Run Variance of the Mean}

\description{
  Convenience function for computing the long-run variance (matrix) of a
  (possibly multivariate) series of observations.
}

\usage{
lrvar(x, type = c("Andrews", "Newey-West"), prewhite = TRUE, adjust = TRUE, \dots)
}

\arguments{
  \item{x}{numeric vector, matrix, or time series.}
  \item{type}{character specifying the type of estimator, i.e., whether
    \code{\link{kernHAC}} for the Andrews quadratic spectral kernel HAC estimator
    is used or \code{\link{NeweyWest}} for the Newey-West Bartlett HAC
    estimator.}
 \item{prewhite}{logical or integer. Should the series be prewhitened?
   Passed to \code{\link{kernHAC}} or \code{\link{NeweyWest}}.}
 \item{adjust}{logical.  Should a finite sample adjustment be made?
   Passed to \code{\link{kernHAC}} or \code{\link{NeweyWest}}.}
 \item{\dots}{further arguments passed on to \code{\link{kernHAC}} or
   \code{\link{NeweyWest}}.}
}

\details{
  \code{lrvar} is a simple wrapper function for computing the long-run variance
  (matrix) of a (possibly multivariate) series \code{x}. First, this simply fits
  a linear regression model \code{x ~ 1} by \code{\link[stats]{lm}}. Second,
  the corresponding variance of the mean(s) is estimated either by \code{\link{kernHAC}}
  (Andrews quadratic spectral kernel HAC estimator) or by \code{\link{NeweyWest}}
  (Newey-West Bartlett HAC estimator).
}

\value{For a univariate series \code{x} a scalar variance is computed. For a
  multivariate series \code{x} the covariance matrix is computed.}

\seealso{\code{\link{kernHAC}}, \code{\link{NeweyWest}}, \code{\link{vcovHAC}}}

\examples{
set.seed(1)
## iid series (with variance of mean 1/n)
## and Andrews kernel HAC (with prewhitening)
x <- rnorm(100)
lrvar(x)

## analogous multivariate case with Newey-West estimator (without prewhitening)
y <- matrix(rnorm(200), ncol = 2)
lrvar(y, type = "Newey-West", prewhite = FALSE)

## AR(1) series with autocorrelation 0.9
z <- filter(rnorm(100), 0.9, method = "recursive")
lrvar(z)
}

\keyword{regression}
\keyword{ts}
