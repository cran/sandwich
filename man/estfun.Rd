\name{estfun}
\alias{estfun}
\alias{estfun.lm}
\alias{estfun.glm}
\alias{estfun.rlm}
\alias{estfun.survreg}
\alias{estfun.coxph}
\alias{estfun.nls}
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
as in \code{\link{coef}} or \code{\link{terms}}, respectively.}

\seealso{\code{\link{lm}}, \code{\link{glm}}}

\references{
  Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \url{http://http://www.jstatsoft.org/v16/i09/}.
}

\examples{
x <- sin(1:10)
y <- rnorm(10)
fm <- lm(y ~ x)

estfun(fm)
residuals(fm) * cbind(1, x)
}

\keyword{regression}
