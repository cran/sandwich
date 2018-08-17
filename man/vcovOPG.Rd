\name{vcovOPG}
\alias{vcovOPG}
\title{Outer-Product-of-Gradients Covariance Matrix Estimation}

\description{
 Outer product of gradients estimation for the covariance matrix of the coefficient estimates
 in regression models.
}

\usage{
vcovOPG(x, adjust = FALSE, \dots)
}

\arguments{
 \item{x}{a fitted model object.}
 \item{adjust}{logical. Should a finite sample adjustment be made?
    This amounts to multiplication with \eqn{n/(n-k)} where \eqn{n} is the
    number of observations and \eqn{k} the number of estimated parameters.}
 \item{\dots}{arguments passed to the \code{\link{estfun}} function.}
}

\details{
In correctly specified models, the \dQuote{meat} matrix (cross product of estimating
functions, see \code{\link{meat}}) and the inverse of the \dQuote{bread} matrix
(inverse of the derivative of the estimating functions, see \code{\link{bread}}) are equal
and correspond to the Fisher information matrix.
Typically, an empirical version of the bread is used for estimation of the information
but alternatively it is also possible to use the meat. This method is also known as
the outer product of gradients (OPG) estimator (Cameron & Trivedi 2005).

Using the \pkg{sandwich} infrastructure, the OPG estimator could easily be computed via
\code{solve(meat(obj))} (modulo scaling). To employ numerically more stable implementation
of the inversion, this simple convenience function can be used: \code{vcovOPG(obj)}.

Note that this only works if the \code{estfun()} method computes the maximum
likelihood scores (and not a scaled version such as least squares scores for
\code{"lm"} objects).
}

\value{
  A matrix containing the covariance matrix estimate.
}

\references{
  Cameron AC and Trivedi PK (2005), \emph{Microeconometrics: Methods and Applications}.
  Cambridge University Press, Cambridge.

  Zeileis A (2006), Object-Oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \url{http://www.jstatsoft.org/v16/i09/}.
}

\seealso{\code{\link{meat}}, \code{\link{bread}}, \code{\link{sandwich}}}

\examples{
## generate poisson regression relationship
x <- sin(1:100)
y <- rpois(100, exp(1 + x))
## compute usual covariance matrix of coefficient estimates
fm <- glm(y ~ x, family = poisson)
vcov(fm)
vcovOPG(fm)
}

\keyword{regression}
\keyword{ts}
