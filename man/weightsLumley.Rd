\name{weightsLumley}
\alias{weightsLumley}
\alias{weave}
\title{Weighted Empirical Adaptive Variance Estimation}
\description{
 A set of functions implementing a class of kernel-based heteroskedasticity
 and autocorrelation consistent (HAC) covariance matrix estimators
 as introduced by Andrews (1991).
}

\usage{
weave(x, order.by = NULL, prewhite = FALSE, C = NULL,
  method = c("truncate", "smooth"), acf = isoacf, adjust = FALSE,
  diagnostics = FALSE, sandwich = TRUE, tol = 1e-7, data = list(), \dots)

weightsLumley(x, order.by = NULL, C = NULL,
  method = c("truncate", "smooth"), acf = isoacf, tol = 1e-7, data = list(), \dots)
}

\arguments{
 \item{x}{a fitted model object.}
 \item{order.by}{Either a vector \code{z} or a formula with a single explanatory
    variable like \code{~ z}. The observations in the model
    are ordered by the size of \code{z}. If set to \code{NULL} (the
    default) the observations are assumed to be ordered (e.g., a
    time series).}
 \item{prewhite}{logical or integer. Should the estimating functions
    be prewhitened? If \code{TRUE} or greater than 0 a VAR model of
    order \code{as.integer(prewhite)} is fitted via \code{ar} with
    method \code{"ols"} and \code{demean = FALSE}.}
 \item{C}{numeric. The cutoff constant \code{C} is by default
    4 for method \code{"truncate"} and 1 for method \code{"smooth"}.}
 \item{method}{a character specifying the method used, see details.}
 \item{acf}{a function that computes the autocorrelation function of 
    a vector, by default \code{\link{isoacf}} is used.}
 \item{adjust}{logical. Should a finite sample adjustment be made?
    This amounts to multiplication with \eqn{n/(n-k)} where \eqn{n} is the
    number of observations and \eqn{k} the number of estimated parameters.}
 \item{diagnostics}{logical. Should additional model diagnostics be returned?
   See \code{\link{vcovHAC}} for details.}
 \item{sandwich}{logical. Should the sandwich estimator be computed?
    If set to \code{FALSE} only the middle matrix is returned.}
 \item{tol}{numeric. Weights that exceed \code{tol} are used for computing
   the covariance matrix, all other weights are treated as 0.}
 \item{data}{an optional data frame containing the variables in the \code{order.by} 
    model. By default the variables are taken from the environment which
    the function is called from.}
 \item{\dots}{currently not used.}
}

\details{
\code{weave} is a convenience interface to \code{\link{vcovHAC}} using 
\code{weightsLumley}: first a weights function is defined and then \code{vcovHAC}
is called.

Both weighting methods are based on some estimate of the autocorrelation
function \eqn{\rho}{r} (as computed by \code{acf}) of the residuals of
the model \code{x}. The weights for the \code{"truncate"} method are 

\deqn{I\{n \rho^2 > C\}}{I{n * r ** 2 > C}}

and the weights for the \code{"smooth"} method are

\deqn{\min\{1, C n \rho^2\}}{min{1, C * n * r ** 2}}

where n is the number of observations in the model an C is the truncation 
constant \code{C}.

Further details can be found in Lumley & Heagerty (1999).
}

\value{
\code{weave} returns the same type of object as \code{\link{vcovHAC}}
which is typically just the covariance matrix.

\code{weightsLumley} returns a vector of weights.
}


\references{
  Lumley A & Heagerty P (1999),
  Weighted Empirical Adaptive Variance Estimators for Correlated Data Regression.
  \emph{Journal of the Royal Statistical Society B}, \bold{61},
  459--477.
}

\seealso{\code{\link{vcovHAC}}, \code{\link{weightsAndrews}},
         \code{\link{kernHAC}}}

\examples{
x <- sin(1:100)
y <- 1 + x + rnorm(100)
fm <- lm(y ~ x)
weave(fm)
vcov(fm)
}

\keyword{regression}
\keyword{ts}
