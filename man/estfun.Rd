\name{estfun}
\alias{estfun}
\alias{estfun.lm}
\alias{estfun.glm}
\alias{estfun.rlm}
\title{Extract Estimating Functions}
\description{
 Generic function for extracting the estimating functions
 of a fitted model with methods for (generalized) linear models.
}
\usage{
estfun(x, \dots)
}
\arguments{
 \item{x}{a fitted model object of class \code{"lm"}, \code{"glm"}
   or \code{"rlm"}.}
 \item{\dots}{currently not used.}
}

\value{A matrix containing the estimating functions.}

\seealso{\code{\link{lm}}, \code{\link{glm}}, \code{\link[MASS]{rlm}}}

\examples{
x <- sin(1:10)
y <- rnorm(10)
fm <- lm(y ~ x)

estfun(fm)
residuals(fm) * cbind(1, x)
}

\keyword{regression}
