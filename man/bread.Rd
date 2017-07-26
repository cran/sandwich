\name{bread}
\alias{bread}
\alias{bread.default}
\alias{bread.lm}
\alias{bread.mlm}
\alias{bread.survreg}
\alias{bread.coxph}
\alias{bread.gam}
\alias{bread.nls}
\alias{bread.rlm}
\alias{bread.hurdle}
\alias{bread.zeroinfl}
\alias{bread.mlogit}
\alias{bread.polr}
\alias{bread.clm}

\title{Bread for Sandwiches}

\description{
 Generic function for extracting an estimator for the bread of
 sandwiches.
}

\usage{
bread(x, \dots)
}

\arguments{
 \item{x}{a fitted model object.}
 \item{\dots}{arguments passed to methods.}
}

\value{
A matrix containing an estimator for the expectation of the negative
derivative of the estimating functions, usually the Hessian.
Typically, this should be an \eqn{k \times k}{k x k} matrix corresponding
to \eqn{k} parameters. The rows and columns should be named
as in \code{\link{coef}} or \code{\link{terms}}, respectively.

The default method tries to extract \code{\link{vcov}} and \code{\link{nobs}}
and simply computes their product.
}

\seealso{\code{\link{lm}}, \code{\link{glm}}}

\references{
  Zeileis A (2006), Object-Oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \url{http://www.jstatsoft.org/v16/i09/}.
}

\examples{
## linear regression
x <- sin(1:10)
y <- rnorm(10)
fm <- lm(y ~ x)

## bread: n * (x'x)^{-1}
bread(fm)
solve(crossprod(cbind(1, x))) * 10
}

\keyword{regression}
