\name{bread}
\alias{bread}
\alias{bread.lm}
\alias{bread.survreg}
\alias{bread.coxph}
\alias{bread.gam}
\alias{bread.nls}
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

\value{A matrix containing an estimator for the expectation of the negative
derivative of the estimating functions.
Typically, this should be an \eqn{k \times k}{k x k} matrix corresponding
to \eqn{k} parameters. The rows and columns should be named
as in \code{\link{coef}} or \code{\link{terms}}, respectively.}

\seealso{\code{\link{lm}}, \code{\link{glm}}}

\references{
  Zeileis A (2006),
  Object-oriented Computation of Sandwich Estimators.
  Package vignette.
}

\examples{
x <- sin(1:10)
y <- rnorm(10)
fm <- lm(y ~ x)

bread(fm)
solve(crossprod(cbind(1, x))) * 10
}

\keyword{regression}
