\name{vcovHAC}
\alias{vcovHAC}
\title{Heteroskedasticity and Autocorrelation Consistent (HAC) Covariance Matrix Estimation}
\description{
 Heteroskedasticity and autocorrelation consistent (HAC) estimation
 of the covariance matrix of the coefficient estimates in a (generalized)
 linear regression model.
}
\usage{
vcovHAC(x, order.by = NULL, prewhite = FALSE, weights = weightsAndrews,
  adjust = TRUE, diagnostics = FALSE, sandwich = TRUE, ar.method = "ols",
  data = list())
}
\arguments{
 \item{x}{a fitted model object of class \code{"lm"} or \code{"glm"}.}
 \item{order.by}{Either a vector \code{z} or a formula with a single explanatory
    variable like \code{~ z}. The observations in the model
    are ordered by the size of \code{z}. If set to \code{NULL} (the
    default) the observations are assumed to be ordered (e.g., a
    time series).}
 \item{prewhite}{logical or integer. Should the estimating functions
    be prewhitened? If \code{TRUE} or greater than 0 a VAR model of
    order \code{as.integer(prewhite)} is fitted via \code{ar} with
    method \code{"ols"} and \code{demean = FALSE}.}
 \item{weights}{Either a vector of weights for the autocovariances or a
    function to compute these weights based on \code{x}, \code{order.by},
    \code{prewhite}, \code{ar.method} and \code{data}. If \code{weights}
    is a function it has to take these arguments. See also details.}
 \item{adjust}{logical. Should a finite sample adjustment be made?
    This amounts to multiplication with $n/(n-k)$ where $n$ is the
    number of observations and $k$ the number of estimated parameters.}
 \item{diagnostics}{logical. Should additional model diagnostics be returned?
    See below for details.}
 \item{sandwich}{logical. Should the sandwich estimator be computed?
    If set to \code{FALSE} only the middle matrix is returned.}
 \item{ar.method}{character. The \code{method} argument passed to
   \code{\link{ar}} for prewhitening.}
 \item{data}{an optional data frame containing the variables in the \code{order.by} 
    model. By default the variables are taken from the environment which
    \code{vcovHAC} is called from.}
}

\details{
This function constructs weighted information sandwich variance estimators
for (generalized) linear models fitted to time series data. These are basically
constructed from weighted sums of autocovariances of the estimation functions
(as extracted by \code{\link{estfun}}). The crucial step is the specification
of weights: the user can either supply \code{vcovHAC} with some vector of 
weights or with a function that computes these weights adaptively (based on
the arguments \code{x}, \code{order.by}, \code{prewhite} and \code{data}). 
Two functions for adaptively choosing weights are implemented in
\code{\link{weightsAndrews}} implementing the results of Andrews (1991) and
in \code{\link{weightsLumley}} implementing the results of Lumley (1999).
The functions \code{\link{kernHAC}} and \code{\link{weave}} respectively
are to more convenient interfaces for \code{vcovHAC} with these functions.

Prewhitening based on VAR approximations is described as suggested in
Andrews & Monahan (1992).

The covariance matrix estimators have been improved by the addition of a bias correction and an
approximate denominator degrees of freedom for test and confidence interval
construction.
}

\value{A matrix containing the covariance matrix estimate. If \code{diagnostics}
was set to \code{TRUE} this has an attribute \code{"diagnostics} which is a list 
with
  \item{bias.correction}{multiplicative bias correction}
  \item{df}{Approximate denominator degrees of freedom}
}


\references{
  Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.

  Andrews DWK & Monahan JC (1992),
  An Improved Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimatior.
  \emph{Econometrica}, \bold{60},
  953--966.

  Lumley A & Heagerty P (1999),
  Weighted Empirical Adaptive Variance Estimators for Correlated Data Regression.
  \emph{Journal of the Royal Statistical Society B}, \bold{61},
  459--477.

  Newey WK & West KD (1987),
  A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix.
  \emph{Econometrica}, \bold{55},
  703--708.
}

\seealso{\code{\link{weightsLumley}}, \code{\link{weightsAndrews}},
         \code{\link{weave}}, \code{\link{kernHAC}}}

\examples{
x <- sin(1:100)
y <- 1 + x + rnorm(100)
fm <- lm(y ~ x)
vcovHAC(fm)
vcov(fm)
}

\keyword{regression}
\keyword{ts}
