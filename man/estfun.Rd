\name{estfun}
\alias{estfun}
\alias{estfun.lm}
\alias{estfun.glm}
\alias{estfun.mlm}
\alias{estfun.rlm}
\alias{estfun.polr}
\alias{estfun.clm}
\alias{estfun.survreg}
\alias{estfun.coxph}
\alias{estfun.nls}
\alias{estfun.hurdle}
\alias{estfun.zeroinfl}
\alias{estfun.mlogit}
\encoding{UTF-8}

\title{Extract Empirical Estimating Functions}

\description{
 Generic function for extracting the empirical estimating functions
 of a fitted model.
}

\usage{
estfun(x, \dots)
}
\arguments{
 \item{x}{a fitted model object.}
 \item{\dots}{arguments passed to methods.}
}

\value{A matrix containing the empirical estimating functions.
Typically, this should be an \eqn{n \times k}{n x k} matrix corresponding
to \eqn{n} observations and \eqn{k} parameters. The columns should be named
as in \code{\link{coef}} or \code{\link{terms}}, respectively.

The estimating function (or score function) for a model is the derivative of the objective function
with respect to the parameter vector. The empirical estimating functions is
the evaluation of the estimating function at the observed data (\eqn{n} observations)
and the estimated parameters (of dimension \eqn{k}).
}

\seealso{\code{\link{lm}}, \code{\link{glm}}}

\references{
Zeileis A (2006).
  \dQuote{Object-Oriented Computation of Sandwich Estimators.}
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  \doi{10.18637/jss.v016.i09}

Zeileis A, Köll S, Graham N (2020).
  \dQuote{Various Versatile Variances: An Object-Oriented Implementation of Clustered Covariances in R.}
  \emph{Journal of Statistical Software}, \bold{95}(1), 1--36.
  \doi{10.18637/jss.v095.i01}
}

\examples{
## linear regression
x <- 1:9
y <- sin(1:9/5)
m <- lm(y ~ x)

## estimating function: (y - x'beta) * x
estfun(m)
residuals(m) * cbind(1, x)
}

\keyword{regression}
