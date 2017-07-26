\name{vcovBS}
\alias{vcovBS}
\alias{vcovBS.default}
\alias{.vcovBSenv}

\title{(Clustered) Bootstrap Covariance Matrix Estimation}

\description{
  Estimation of basic bootstrap covariances, using a simple
  (clustered) case-based resampling.
  covariance matrices using an object-oriented approach.
}

\usage{
vcovBS(x, \dots)

\method{vcovBS}{default}(x, cluster = NULL, R = 250, start = FALSE, \dots)
}

\arguments{
  \item{x}{a fitted model object.}
  \item{cluster}{a variable indicating the clustering of observations.
    By default each observation is its own cluster.}
  \item{R}{integer. Number of bootstrap replications.}
  \item{start}{logical. Should \code{coef(x)} be passed as \code{start}
    to the \code{update(x, subset = ...)} call? In case the model \code{x}
    is computed by some numeric iteration, this may speed up the bootstrapping.}
  \item{\dots}{arguments passed to methods (to \code{update}
    in case of the default method).}
}

\details{
  Basic (clustered) bootstrap covariance matrix estimation is provided by
  \code{vcovBS}. The default method samples clusters (where each observation
  is its own cluster by default), i.e. case-based resampling. For obtaining
  a covariance matrix estimate it is assumed that an \code{\link[stats]{update}}
  of the model with the resampled subset can be obtained, the \code{\link[stats]{coef}}
  extracted, and finally the covariance computed with \code{\link[stats]{cov}}.
  
  Note: The \code{update} model is evaluated in the \code{environment(terms(x))}.
}

\value{
  A matrix containing the covariance matrix estimate.
}

\seealso{\code{\link{vcovCL}}}

\examples{
## Petersen's data
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

## comparison of different standard errors
set.seed(1)
cbind(
  "classical" = sqrt(diag(vcov(m))),
  "HC-cluster" = sqrt(diag(vcovCL(m, cluster = PetersenCL$firm))),
  "BS-cluster" = sqrt(diag(vcovBS(m, cluster = PetersenCL$firm)))
)
}

\keyword{regression}
