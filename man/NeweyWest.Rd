\name{NeweyWest}
\alias{bwNeweyWest}
\alias{NeweyWest}
\title{Newey-West HAC Covariance Matrix Estimation}
\description{
 A set of functions implementing the Newey & West (1987, 1994) heteroskedasticity
 and autocorrelation consistent (HAC) covariance matrix estimators.
}

\usage{
NeweyWest(x, lag = NULL, order.by = NULL, prewhite = TRUE, adjust = FALSE, 
  diagnostics = FALSE, sandwich = TRUE, ar.method = "ols", data = list(),
  verbose = FALSE)

bwNeweyWest(x, order.by = NULL, kernel = c("Bartlett", "Parzen",
  "Quadratic Spectral", "Truncated", "Tukey-Hanning"), weights = NULL,
  prewhite = 1, ar.method = "ols", data = list(), \dots)
}

\arguments{
 \item{x}{a fitted model object of class \code{"lm"} or \code{"glm"}.}
 \item{lag}{integer specifying the maximum lag with positive 
    weight for the Newey-West estimator. If set to \code{NULL}
    \code{floor(bwNeweyWest(x, ...))} is used.}
 \item{order.by}{Either a vector \code{z} or a formula with a single explanatory
    variable like \code{~ z}. The observations in the model
    are ordered by the size of \code{z}. If set to \code{NULL} (the
    default) the observations are assumed to be ordered (e.g., a
    time series).}
 \item{prewhite}{logical or integer. Should the estimating functions
    be prewhitened? If \code{TRUE} or greater than 0 a VAR model of
    order \code{as.integer(prewhite)} is fitted via \code{ar} with
    method \code{"ols"} and \code{demean = FALSE}. The default is to
    use VAR(1) prewhitening.}
 \item{kernel}{a character specifying the kernel used. All kernels used
    are described in Andrews (1991). \code{bwNeweyWest} can only
    compute bandwidths for \code{"Bartlett"}, \code{"Parzen"} and
    \code{"Quadratic Spectral"}.}
 \item{adjust}{logical. Should a finite sample adjustment be made?
    This amounts to multiplication with $n/(n-k)$ where $n$ is the
    number of observations and $k$ the number of estimated parameters.}
 \item{diagnostics}{logical. Should additional model diagnostics be returned?
   See \code{\link{vcovHAC}} for details.}
 \item{sandwich}{logical. Should the sandwich estimator be computed?
    If set to \code{FALSE} only the middle matrix is returned.}
 \item{ar.method}{character. The \code{method} argument passed to
   \code{\link{ar}} for prewhitening (only, not for bandwidth selection).}
 \item{data}{an optional data frame containing the variables in the \code{order.by} 
    model. By default the variables are taken from the environment which
    the function is called from.}
 \item{verbose}{logical. Should the lag truncation parameter used be
    printed?}
 \item{weights}{numeric. A vector of weights used for weighting the estimated
   coefficients of the approximation model (as specified by \code{approx}). By
   default all weights are 1 except that for the intercept term (if there is more than
   one variable).}
 \item{\dots}{currently not used.}
}

\details{\code{NeweyWest} is a convenience interface to \code{\link{vcovHAC}} using 
Bartlett kernel weights as described in Newey & West (1987, 1994).
The automatic bandwidth selection procedure described in Newey & West (1994)
is used as the default and can also be supplied to \code{kernHAC} for the
Parzen and quadratic spectral kernel. It is implemented in \code{bwNeweyWest}
which does not truncate its results - if the results for the Parzen and Bartlett
kernels should be truncated, this has to be applied afterwards. For Bartlett 
weights this is implemented in \code{NeweyWest}.

To obtain the estimator described in Newey & West (1987), prewhitening has to
be suppressed.
}

\value{
\code{NeweyWest} returns the same type of object as \code{\link{vcovHAC}}
which is typically just the covariance matrix.

\code{bwNeweyWest} returns the selected bandwidth parameter.
}


\references{
  Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.

  Newey WK & West KD (1987),
  A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix.
  \emph{Econometrica}, \bold{55},
  703--708.

  Newey WK & West KD (1994),
  Automatic Lag Selection in Covariance Matrix Estimation.
  \emph{Review of Economic Studies}, \bold{61},
  631--653.
}

\seealso{\code{\link{vcovHAC}}, \code{\link{weightsAndrews}},
         \code{\link{kernHAC}}}

\examples{
## fit investment equation
data(Investment)
fm <- lm(RealInv ~ RealGNP + RealInt, data = Investment)

## Newey & West (1994) compute this type of estimator
NeweyWest(fm)

## The Newey & West (1987) estimator requires specification
## of the lag and suppression of prewhitening
NeweyWest(fm, lag = 4, prewhite = FALSE)

## bwNeweyWest() can also be passed to kernHAC(), e.g.
## for the quadratic spectral kernel
kernHAC(fm, bw = bwNeweyWest)
}

\keyword{regression}
\keyword{ts}
