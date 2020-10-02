\name{meat}
\alias{meat}
\encoding{UTF-8}

\title{A Simple Meat Matrix Estimator}

\description{
 Estimating the variance of the estimating functions of
 a regression model by cross products of the empirical
 estimating functions.
}
 
\usage{
meat(x, adjust = FALSE, \dots)
}
\arguments{
 \item{x}{a fitted model object.}
 \item{adjust}{logical. Should a finite sample adjustment be made?
    This amounts to multiplication with \eqn{n/(n-k)} where \eqn{n} is the
    number of observations and \eqn{k} the number of estimated parameters.}
 \item{\dots}{arguments passed to the \code{\link{estfun}} function.}
}

\details{For some theoretical background along with implementation
details see Zeileis (2006).
}

\value{A \eqn{k \times k}{k x k} matrix corresponding containing
the scaled cross products of the empirical estimating functions.}

\seealso{\code{\link{sandwich}}, \code{\link{bread}}, \code{\link{estfun}}}

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
x <- sin(1:10)
y <- rnorm(10)
fm <- lm(y ~ x)

meat(fm)
meatHC(fm, type = "HC")
meatHAC(fm)
}

\keyword{regression}
