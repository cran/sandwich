\name{vcovBS}
\alias{vcovBS}
\alias{vcovBS.default}
\alias{vcovBS.lm}
\alias{vcovBS.glm}
\alias{.vcovBSenv}
\encoding{UTF-8}

\title{(Clustered) Bootstrap Covariance Matrix Estimation}

\description{
  Object-oriented estimation of basic bootstrap covariances, using
  simple (clustered) case-based resampling, plus more refined methods
  for \code{\link[stats]{lm}} and \code{\link[stats]{glm}} models.
}

\usage{
vcovBS(x, \dots)

\method{vcovBS}{default}(x, cluster = NULL, R = 250, start = FALSE, \dots,
  fix = FALSE, use = "pairwise.complete.obs", applyfun = NULL, cores = NULL)

\method{vcovBS}{lm}(x, cluster = NULL, R = 250, type = "xy", \dots,
  fix = FALSE, use = "pairwise.complete.obs", applyfun = NULL, cores = NULL)

\method{vcovBS}{glm}(x, cluster = NULL, R = 250, start = FALSE, \dots,
  fix = FALSE, use = "pairwise.complete.obs", applyfun = NULL, cores = NULL)
}

\arguments{
  \item{x}{a fitted model object.}
  \item{cluster}{a variable indicating the clustering of observations,
    a \code{list} (or \code{data.frame}) thereof, or a formula specifying
    which variables from the fitted model should be used (see examples).
    By default (\code{cluster = NULL}), either \code{attr(x, "cluster")} is used
    (if any) or otherwise every observation is assumed to be its own cluster.}
  \item{R}{integer. Number of bootstrap replications.}
  \item{start}{logical. Should \code{coef(x)} be passed as \code{start}
    to the \code{update(x, subset = ...)} call? In case the model \code{x}
    is computed by some numeric iteration, this may speed up the bootstrapping.}
  \item{type}{character or function. The character string specifies the type of
    bootstrap to use: One of \code{"xy"}, \code{"residual"}, \code{"wild"}
    (or equivalently: \code{"wild-rademacher"} or \code{"rademacher"}),
    \code{"mammen"} (or \code{"wild-mammen"}), \code{"norm"}
    (or \code{"wild-norm"}), \code{"webb"} (or \code{"wild-webb"}).
    Alternatively, \code{type} can be a \code{function(n)} for drawing
    wild bootstrap factors.}
  \item{\dots}{arguments passed to methods. For the default method, this is
    passed to \code{update}, and for the \code{lm} method to \code{lm.fit}.}
  \item{fix}{logical. Should the covariance matrix be fixed to be
    positive semi-definite in case it is not?}
  \item{use}{character. Specification passed to \code{\link[stats]{cov}} for
    handling missing coefficients/parameters.}
  \item{applyfun}{an optional \code{\link[base]{lapply}}-style function with arguments
    \code{function(X, FUN, \dots)}. It is used for refitting the model to the
    bootstrap samples. The default is to use the basic \code{lapply}
    function unless the \code{cores} argument is specified (see below).}
  \item{cores}{numeric. If set to an integer the \code{applyfun} is set to    
    \code{\link[parallel]{mclapply}} with the desired number of \code{cores},
    except on Windows where \code{\link[parallel]{parLapply}} with
    \code{makeCluster(cores)} is used.}
}

\details{
  Clustered sandwich estimators are used to adjust inference when errors
  are correlated within (but not between) clusters. See the documentation for \code{\link{vcovCL}} 
  for specifics about covariance clustering. This function allows
  for clustering in arbitrary many cluster dimensions (e.g., firm, time, industry), given all
  dimensions have enough clusters (for more details, see Cameron et al. 2011).
  Unlike \code{vcovCL}, \code{vcovBS} uses a bootstrap rather than an asymptotic solution.

  Basic (clustered) bootstrap covariance matrix estimation is provided by
  the default \code{vcovBS} method. It samples clusters (where each observation
  is its own cluster by default), i.e., using case-based resampling. For obtaining
  a covariance matrix estimate it is assumed that an \code{\link[stats]{update}}
  of the model with the resampled \code{subset} can be obtained, the \code{\link[stats]{coef}}
  extracted, and finally the covariance computed with \code{\link[stats]{cov}}.
  
  The \code{update} model is evaluated in the \code{environment(terms(x))} (if available).
  To speed up computations two further arguments can be leveraged.
  \enumerate{
    \item Instead of \code{\link[base]{lapply}} a parallelized function such as
      \code{\link[parallel]{parLapply}} or \code{\link[parallel]{mclapply}}
      can be specified to iterate over the bootstrap replications. For the latter,
      specifying \code{cores = ...} is a convenience shortcut.
    \item When specifying \code{start = TRUE}, the \code{coef(x)} are passed to
      \code{update} as \code{start = coef(x)}. This may not be supported by all
      model fitting functions and is hence not turned on by default.
  }
  
  The ``xy'' or ``pairs'' bootstrap is consistent for heteroscedasticity and clustered errors, 
  and converges to the asymptotic solution used in \code{vcovCL}
  as \code{R}, \eqn{n}, and \eqn{g} become large (\eqn{n} and \eqn{g} are the number of
  observations and the number of clusters, respectively; see Efron 1979, or Mammen 1992, for a 
  discussion of bootstrap asymptotics). For small \eqn{g}--particularly under 30 groups--the
  bootstrap will converge to a slightly different value than the asymptotic method, due to
  the limited number of distinct bootstrap replications possible (see Webb 2014 for a discussion
  of this phenomonon). The bootstrap will not necessarily converge to an asymptotic estimate
  that has been corrected for small samples.
  
  The xy approach to bootstrapping is generally only of interest to the
  practitioner when the asymptotic solution is unavailable (this can happen when using
  estimators that have no \code{estfun} function, for example). The residual bootstrap,
  by contrast, is rarely of practical interest, because while it provides consistent
  inference for for clustered standard errors, it is not robust to heteroscedasticity.
  More generally, bootstrapping is useful when the bootstrap makes different assumptions than the asymptotic
  estimator, in particular when the number of clusters is small and large \eqn{n} or
  \eqn{g} assumptions are unreasonable. Bootstrapping is also often effective for nonlinear models,
  particularly in smaller samples, where asymptotic approaches often perform relatively poorly.
  See Cameron and Miller (2015) for further discussion of bootstrap techniques in practical applications,
  and Zeileis et al. (2020) show simulations comparing \code{vcovBS} to \code{vcovCL} in several
  settings.

  The \code{\link[stats]{glm}} method works essentially like the default method but call
  \code{\link[stats]{glm.fit}} instead of code{update}.
  
  The \code{\link[stats]{lm}} method provides additional bootstrapping \code{type}s
  and computes the bootstrapped coefficient estimates somewhat more efficiently using
  \code{\link[stats]{lm.fit}} (for case-based resampling) or \code{\link[base]{qr.coef}}
  rather than \code{update}. The default \code{type} is case-based resampling
  (\code{type = "xy"}) as in the default method. Alternative \code{type} specifications are:
  \itemize{
    \item \code{"residual"}. The residual cluster bootstrap resamples the residuals (as above,
      by cluster) which are subsequently added to the fitted values to obtain the bootstrapped
      response variable: \eqn{y^{*} = \hat{y} + e^{*}}{y* = yhat + e*}.
      Coefficients can then be estimated using \code{qr.coef()}, reusing the
      QR decomposition from the original fit. As Cameron et al. (2008) point out,
      the residual cluster bootstrap is not well-defined when the clusters are unbalanced as 
      residuals from one cluster cannot be easily assigned to another cluster with different size.
      Hence a warning is issued in that case.
    \item \code{"wild"} (or equivalently \code{"wild-rademacher"} or \code{"rademacher"}).
      The wild cluster bootstrap does not actually resample the residuals but instead reforms the
      dependent variable by multiplying the residual by a randomly drawn value and adding the
      result to the fitted value: \eqn{y^{*} = \hat{y} + e \cdot w}{y* = yhat + e * w}
      (see Cameron et al. 2008). By default, the factors are drawn from the Rademacher distribution:
      \code{function(n) sample(c(-1, 1), n, replace = TRUE)}.
    \item \code{"mammen"} (or \code{"wild-mammen"}). This draws the wild bootstrap factors as
      suggested by Mammen (1993):
      \code{sample(c(-1, 1) * (sqrt(5) + c(-1, 1))/2, n, replace = TRUE, prob = (sqrt(5) + c(1, -1))/(2 * sqrt(5)))}.
    \item \code{"webb"} (or \code{"wild-webb"}). This implements the six-point distribution
      suggested by Webb (2014), which may improve inference when the number of clusters is small:
      \code{sample(c(-sqrt((3:1)/2), sqrt((1:3)/2)), n, replace = TRUE)}.
    \item \code{"norm"} (or \code{"wild-norm"}). The standard normal/Gaussian distribution
      is used for drawing the wild bootstrap factors: \code{function(n) rnorm(n)}.
    \item User-defined function. This needs of the form as above, i.e., a \code{function(n)}
      returning a vector of random wild bootstrap factors of corresponding length.
   }
}

\value{
  A matrix containing the covariance matrix estimate.
}

\references{
Cameron AC, Gelbach JB, Miller DL (2008).
  \dQuote{Bootstrap-Based Improvements for Inference with Clustered Errors},
  \emph{The Review of Economics and Statistics}, \bold{90}(3), 414--427.
  \doi{10.3386/t0344}
  
Cameron AC, Gelbach JB, Miller DL (2011).
  \dQuote{Robust Inference with Multiway Clustering},
  \emph{Journal of Business & Economic Statistics}, \bold{29}(2), 238--249.
  \doi{10.1198/jbes.2010.07136}
  
Cameron AC, Miller DL (2015). 
  \dQuote{A Practitioner's Guide to Cluster-Robust Inference}, 
  \emph{Journal of Human Resources}, \bold{50}(2), 317--372.
  \doi{10.3368/jhr.50.2.317}
  
Efron B (1979). 
  \dQuote{Bootstrap Methods: Another Look at the Jackknife}, 
  \emph{The Annals of Statistics}, \bold{7}(1), 1--26.
  \doi{10.1214/aos/1176344552}
  
Mammen, E (1992). 
  \dQuote{When Does Bootstrap Work?: Asymptotic Results and Simulations}, 
  \emph{Lecture Notes in Statistics}, \bold{77}.
  Springer Science & Business Media.
  
Mammen, E (1993). 
  \dQuote{Bootstrap and Wild Bootstrap for High Dimensional Linear Models}, 
  \emph{The Annals of Statistics}, \bold{21}(1), 255--285.
  \doi{10.1214/aos/1176349025}
  
Webb, MD (2014).
  \dQuote{Reworking Wild Bootstrap Based Inference for Clustered Errors},
  Working Paper 1315, \emph{Queen's Economics Department.}
  \url{http://qed.econ.queensu.ca/working_papers/papers/qed_wp_1315.pdf}.

Zeileis A, Köll S, Graham N (2020).
  \dQuote{Various Versatile Variances: An Object-Oriented Implementation of Clustered Covariances in R.}
  \emph{Journal of Statistical Software}, \bold{95}(1), 1--36.
  \doi{10.18637/jss.v095.i01}
}

\seealso{\code{\link{vcovCL}}}

\examples{
## Petersen's data
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

## comparison of different standard errors
suppressWarnings(RNGversion("3.5.0"))
set.seed(1)
cbind(
  "classical" = sqrt(diag(vcov(m))),
  "HC-cluster" = sqrt(diag(vcovCL(m, cluster = ~ firm))),
  "BS-cluster" = sqrt(diag(vcovBS(m, cluster = ~ firm)))
)

## two-way wild cluster bootstrap with Mammen distribution
vcovBS(m, cluster = ~ firm + year, type = "wild-mammen")
}

\keyword{regression}
\keyword{bootstrap}
