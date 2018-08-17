\name{vcovHC}
\alias{vcovHC}
\alias{vcovHC.default}
\alias{vcovHC.mlm}
\alias{meatHC}
\title{Heteroskedasticity-Consistent Covariance Matrix Estimation}
\description{
 Heteroskedasticity-consistent estimation of the covariance matrix of the
 coefficient estimates in regression models.
}
\usage{
vcovHC(x, \dots)

\method{vcovHC}{default}(x,
  type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5"),
  omega = NULL, sandwich = TRUE, \dots)

meatHC(x, type = , omega = NULL, \dots)
}
\arguments{
 \item{x}{a fitted model object.}
 \item{type}{a character string specifying the estimation type. For
     details see below.}
 \item{omega}{a vector or a
     function depending on the arguments \code{residuals}
     (the working residuals of the model), \code{diaghat} (the diagonal 
     of the corresponding hat matrix) and \code{df} (the residual degrees of
     freedom). For details see below.}
 \item{sandwich}{logical. Should the sandwich estimator be computed?
    If set to \code{FALSE} only the meat matrix is returned.}
 \item{\dots}{arguments passed to \code{\link{sandwich}} (in \code{vcovHC})
   and \code{\link{estfun}} (in \code{meatHC}), respectively.}
}

\details{The function \code{meatHC} is the real work horse for estimating
  the meat of HC sandwich estimators -- the default \code{vcovHC} method
  is a wrapper calling
  \code{\link{sandwich}} and \code{\link{bread}}. See Zeileis (2006) for
  more implementation details. The theoretical background, exemplified
  for the linear regression model, is described below and in Zeileis (2004).
  Analogous formulas are employed for other types of models.

  When \code{type = "const"} constant variances are assumed and
  and \code{vcovHC} gives the usual estimate of the covariance matrix of
  the coefficient estimates:

  \deqn{\hat \sigma^2 (X^\top X)^{-1}}{sigma^2 (X'X)^{-1}}

  All other methods do not assume constant variances and are suitable in case of
  heteroskedasticity. \code{"HC"} (or equivalently \code{"HC0"}) gives White's
  estimator, the other estimators are refinements of this. They are all of form
  
  \deqn{(X^\top X)^{-1} X^\top \Omega X (X^\top X)^{-1}}{(X'X)^{-1} X' Omega X (X'X)^{-1}}

  and differ in the choice of Omega. This is in all cases a diagonal matrix whose 
  elements can be either supplied as a vector \code{omega} or as a
  a function \code{omega} of the residuals, the diagonal elements of the hat matrix and
  the residual degrees of freedom. For White's estimator
  
  \code{omega <- function(residuals, diaghat, df) residuals^2}
  
  Instead of specifying the diagonal \code{omega} or a function for
  estimating it, the \code{type} argument can be used to specify the 
  HC0 to HC5 estimators. If \code{omega} is used, \code{type} is ignored.
  
  Long \& Ervin (2000) conduct a simulation study of HC estimators (HC0 to HC3) in
  the linear regression model, recommending to use HC3 which is thus the
  default in \code{vcovHC}. Cribari-Neto (2004), Cribari-Neto, Souza, \& Vasconcellos (2007),
  and Cribari-Neto \& Da Silva (2011), respectively, suggest the HC4, HC5, and
  modified HC4m type estimators. All of them are tailored to take into account
  the effect of leverage points in the design matrix. For more details see the references.
}

\value{
  A matrix containing the covariance matrix estimate.
}

\references{
Cribari-Neto F. (2004), Asymptotic Inference under Heteroskedasticity
of Unknown Form. \emph{Computational Statistics & Data Analysis}
\bold{45}, 215--233.

Cribari-Neto F., Da Silva W.B. (2011), A New Heteroskedasticity-Consistent
Covariance Matrix Estimator for the Linear Regression Model.
\emph{Advances in Statistical Analysis}, \bold{95}(2), 129--146.

Cribari-Neto F., Souza T.C., Vasconcellos, K.L.P. (2007), Inference under
Heteroskedasticity and Leveraged Data. \emph{Communications in Statistics -- Theory and
Methods}, \bold{36}, 1877--1888. Errata: \bold{37}, 3329--3330, 2008.

Long J. S., Ervin L. H. (2000), Using Heteroscedasticity Consistent Standard
Errors in the Linear Regression Model. \emph{The American Statistician},
\bold{54}, 217--224.

MacKinnon J. G., White H. (1985), Some Heteroskedasticity-Consistent
Covariance Matrix Estimators with Improved Finite Sample Properties.
\emph{Journal of Econometrics} \bold{29}, 305--325.

White H. (1980), A Heteroskedasticity-Consistent Covariance Matrix and
a Direct Test for Heteroskedasticity. \emph{Econometrica} \bold{48},
817--838.

Zeileis A (2004), Econometric Computing with HC and HAC Covariance Matrix
Estimators. \emph{Journal of Statistical Software}, \bold{11}(10), 1--17.
URL \url{http://www.jstatsoft.org/v11/i10/}.

Zeileis A (2006), Object-Oriented Computation of Sandwich Estimators.
\emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
URL \url{http://www.jstatsoft.org/v16/i09/}.
}

\seealso{\code{\link{lm}}, \code{\link[car]{hccm}},
\code{\link[lmtest]{bptest}}, \code{\link[car]{ncv.test}}}

\examples{
## generate linear regression relationship
## with homoskedastic variances
x <- sin(1:100)
y <- 1 + x + rnorm(100)
## model fit and HC3 covariance
fm <- lm(y ~ x)
vcovHC(fm)
## usual covariance matrix
vcovHC(fm, type = "const")
vcov(fm)

sigma2 <- sum(residuals(lm(y ~ x))^2)/98
sigma2 * solve(crossprod(cbind(1, x)))
}

\keyword{regression}
\keyword{ts}
