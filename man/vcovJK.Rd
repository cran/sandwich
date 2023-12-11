\name{vcovJK}
\alias{vcovJK}
\alias{vcovJK.default}
\encoding{UTF-8}

\title{(Clustered) Jackknife Covariance Matrix Estimation}

\description{
  Object-oriented estimation of jackknife covariances, i.e., based on
  the centered outer product of leave-on-out estimates of the model
  coefficients/parameters.
}

\usage{
vcovJK(x, \dots)

\method{vcovJK}{default}(x, cluster = NULL, center = "mean", \dots)
}

\arguments{
  \item{x}{a fitted model object.}
  \item{cluster}{a variable indicating the clustering of observations,
    a \code{list} (or \code{data.frame}) thereof, or a formula specifying
    which variables from the fitted model should be used (see examples).
    By default (\code{cluster = NULL}), either \code{attr(x, "cluster")} is used
    (if any) or otherwise every observation is assumed to be its own cluster.}
  \item{center}{character specifying how to center the coefficients from
    all jacknife samples (each dropping one observational unit/cluster).
    By default the coefficients are centered by their \code{"mean"} across the
    sample or, alternatively, by the original full-sample \code{"estimate"}.}
  \item{\dots}{arguments passed to methods. For the default method, this is
    passed to \code{\link{vcovBS}}.}
}

\details{
  Jackknife covariance estimation is based on leave-one-out estimates of the
  coefficients/parameters of a model. This means that the model is reestimated
  after dropping each observational unit once, i.e., each individual observation
  in independent observations or each cluster in dependent data. The covariance
  matrix is then constructed from the scaled outer product of the centered
  jackknife estimates. Centering can either be done by the mean of the jackknife
  coefficients (default) or by the original full-sample estimates. Scaling is done
  by (N - 1)/N where N is the number of observational units.
  
  Recent research has shown that the jackknife covariance estimate have particularly
  useful properties in practice: they are not downward biased and yield better
  coverage rates for confidence intervals compared to other "robust" covariance
  estimates. See MacKinnon et al. (2022) and Hansen (2022) for more details.
  
  As jackknife covariances are also based on reestimation of the coefficients on
  subsamples, their computation is very similar to bootstrap covariances. Hence,
  the \code{\link{vcovBS}} methods provided in the package all offer an argument
  \code{vcovBS(..., type = "jackknife")}. This is called by the default
  \code{vcovJK} method. Therefore, see the arguments of \code{vcovBS} for further
  details, e.g., for leveraging multicore computations etc.

  In the linear regression model, the jackknife covariance can actually be computed
  without reestimating the coefficients but using only the full-sample estimates and
  certain elements of the so-called hat matrix. Namly the diagonal elements or
  blocks of elements from the hat matrix are needed for independent observations and
  clustered data, respectively. These alternative computations of the jackknife
  covariances are available in \code{\link{vcovHC}} and \code{\link{vcovCL}}, respectively,
  in both cases with argument \code{type = "HC3"}. To obtain HC3 covariances that exactly
  match the jackknife covariances, the jackknife has to be centered with the full-sample
  estimates and the right finite-sample adjustment has to be selected for the HC3.
  
  In small to moderate sample sizes, the HC3 estimation techniques are typically much
  faster than the jackknife. However, in large samples it may become impossible to
  compute the HC3 covariances while the jackknife approach is still feasible.
}

\value{
  A matrix containing the covariance matrix estimate.
}

\references{
Bell RM, McCaffrey DF (2002).
  \dQuote{Bias Reduction in Standard Errors for Linear Regression with Multi-Stage Samples},
  \emph{Survey Methodology}, \bold{28}(2), 169--181.

Hansen BE (2022).
  \dQuote{Jackknife Standard Errors for Clustered Regression}, Working Paper, August 2022.
  \url{https://www.ssc.wisc.edu/~bhansen/papers/tcauchy.html}

MacKinnon JG, Nielsen MØ, Webb MD (2022). 
  \dQuote{Cluster-Robust Inference: A Guide to Empirical Practice}, 
  \emph{Journal of Econometrics}, Forthcoming.
  \doi{10.1016/j.jeconom.2022.04.001}
  
Zeileis A, Köll S, Graham N (2020).
  \dQuote{Various Versatile Variances: An Object-Oriented Implementation of Clustered Covariances in R.}
  \emph{Journal of Statistical Software}, \bold{95}(1), 1--36.
  \doi{10.18637/jss.v095.i01}
}

\seealso{\code{\link{vcovJK}}, \code{\link{vcovHC}}, \code{\link{vcovCL}}}

\examples{
## cross-section data
data("PublicSchools", package = "sandwich")
m1 <- lm(Expenditure ~ poly(Income, 2), data = PublicSchools)
vcovJK(m1, center = "estimate")
vcovHC(m1, type = "HC3") * (nobs(m1) - 1)/nobs(m1)

## clustered data
data("PetersenCL", package = "sandwich")
m2 <- lm(y ~ x, data = PetersenCL)

## jackknife estimator coincides with HC3 (aka CV3)
vcovJK(m2, cluster = ~ firm, center = "estimate")
vcovCL(m2, cluster = ~ firm, type = "HC3", cadjust = FALSE)
}

\keyword{regression}
\keyword{bootstrap}
