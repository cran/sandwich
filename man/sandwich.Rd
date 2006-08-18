\name{sandwich}
\alias{sandwich}
\title{Making Sandwiches with Bread and Meat}
\description{
 Constructing sandwich covariance matrix estimators by
 multiplying bread and meat matrices.
}
\usage{
sandwich(x, bread. = bread, meat. = meat, \dots)
}
\arguments{
 \item{x}{a fitted model object.}
 \item{bread.}{either a bread matrix or a function for computing
   this via \code{bread.(x)}.}
 \item{meat.}{either a bread matrix or a function for computing
   this via \code{meat.(x, \dots)}.}
 \item{\dots}{arguments passed to the \code{meat} function.}
}

\details{\code{sandwich} is a simple convenience function that
takes a bread matrix (i.e., estimator of the expectation of the negative
derivative of the estimating functions) and a meat matrix (i.e.,
estimator of the variance of the estimating functions) and multiplies
them to a sandwich with meat between two slices of bread. By default
\code{\link{bread}} and \code{\link{meat}} are called. 

Some theoretical background along with implementation details is given
in Zeileis (2006).
}

\value{A matrix containing the sandwich covariance matrix estimate.
Typically, this should be an \eqn{k \times k}{k x k} matrix corresponding
to \eqn{k} parameters.}

\seealso{\code{\link{bread}}, \code{\link{meat}}, \code{\link{meatHC}}, \code{\link{meatHAC}}}

\references{
  Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \url{http://http://www.jstatsoft.org/v16/i09/}.
}

\examples{
x <- sin(1:10)
y <- rnorm(10)
fm <- lm(y ~ x)

sandwich(fm)
vcovHC(fm, type = "HC")
}

\keyword{regression}
