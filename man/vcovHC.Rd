\name{vcovHC}
\alias{vcovHC}
\title{Heteroskedasticity-Consistent Covariance Matrix Estimation}
\description{
 Heteroskedasticity-consistent estimation of the covariance matrix of the
 coefficient estimates in a linear regression model.
}
\usage{
vcovHC(x, order.by = NULL, data = list(),
  type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4"),
  omega = NULL)
}
\arguments{
 \item{x}{a fitted model object of class \code{"lm"}.}
 \item{order.by}{formula. A formula with a single explanatory
    variable like \code{~ x}. The observations in the model
    are ordered by the size of \code{x}. If set to \code{NULL} (the
    default) the observations are assumed to be ordered (e.g. a
    time series).}
 \item{data}{an optional data frame containing the variables in the \code{order.by} 
    model. By default the variables are taken from the environment which
    \code{vcovHC} is called from.}
 \item{type}{a character string specifying the estimation type. For
     details see below.}
 \item{omega}{a function depending on the arguments \code{residuals}
     (the residuals of the linear model), \code{diaghat} (the diagonal 
     of the corresponding hat matrix) and \code{df} (the residual degrees of
     freedom). For details see below.}
}
\details{When \code{type = "const"} constant variances are assumed and
  and \code{covHC} gives the usual estimate of the covariance matrix of
  the coefficient estimates:

  \deqn{\hat \sigma^2 (X^\top X)^{-1}}{sigma^2 (X'X)^{-1}}

  All other methods do not assume constant variances and are suitable in case of
  heteroskedasticity. \code{"HC"} (or equivalently \code{"HC0"}) gives White's
  estimator, the other estimators are refinements of this. They are all of form
  
  \deqn{(X^\top X)^{-1} X^\top \Omega X (X^\top X)^{-1}}{(X'X)^{-1} X' Omega X (X'X)^{-1}}

  and differ in the choice of Omega. This is in all cases a diagonal matrix whose 
  elements are function of the residuals, the diagonal elements of the hat matrix and
  the residual degrees of freedom. For White's estimator
  
  \code{omega <- function(residuals, diaghat, df) residuals^2}
  
  Instead of specifying a \code{type} the argument \code{omega} can also be specified
  directly. For details see the references.
}
\value{A matrix containing the covariance matrix estimate.
}
\references{

Cribari-Neto F. (2004), Asymptotic inference under heteroskedasticity
of unknown form. \emph{Computational Statistics \& Data Analysis}
\bold{45}, 215-233.

MacKinnon J. G., White H. (1985),
Some heteroskedasticity-consistent
covariance matrix estimators with improved finite sample properties.
\emph{Journal of Econometrics} \bold{29}, 305-325.

White H. (1980), A heteroskedasticity-consistent covariance matrix and
a direct test for heteroskedasticity. \emph{Econometrica} \bold{48},
817-838.

}

\seealso{\code{\link{lm}}, \code{\link[car]{hccm}},
\code{\link[lmtest]{bptest}}, \code{\link[car]{ncv.test}}}

\examples{
## generate linear regression relationship
## with homoskedastic variances
x <- sin(1:100)
y <- 1 + x + rnorm(100)
## compute usual covariance matrix of coefficient estimates
fm <- lm(y ~ x)
vcovHC(fm, type="const")
vcov(fm)

sigma2 <- sum(residuals(lm(y~x))^2)/98
sigma2 * solve(crossprod(cbind(1,x)))
}

\keyword{regression}
\keyword{ts}
