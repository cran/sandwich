\name{weightsAndrews}
\alias{weightsAndrews}
\alias{bwAndrews}
\alias{kernHAC}
\title{Kernel-based HAC Covariance Matrix Estimation}
\description{
 A set of functions implementing a class of kernel-based heteroskedasticity
 and autocorrelation consistent (HAC) covariance matrix estimators
 as introduced by Andrews (1991).
}

\usage{
kernHAC(x, order.by = NULL, prewhite = 1, bw = NULL,
  kernel = c("Quadratic Spectral", "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),
  approx = c("AR(1)", "ARMA(1,1)"), diagnostics = FALSE, sandwich = TRUE,
  ar.method = "ols", tol = 1e-7, data = list(), \dots)

weightsAndrews(x, order.by = NULL, bw = NULL,
  kernel = c("Quadratic Spectral", "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),
  prewhite = 1, ar.method = "ols", tol = 1e-7, data = list(), \dots)

bwAndrews(x, order.by = NULL, kernel = c("Quadratic Spectral", "Truncated",
  "Bartlett", "Parzen", "Tukey-Hanning"), approx = c("AR(1)", "ARMA(1,1)"),
  weights = NULL, prewhite = 1, ar.method = "ols", data = list())
}

\arguments{
 \item{x}{a fitted model object of class \code{"lm"} or \code{"glm"}.}
 \item{order.by}{Either a vector \code{z} or a formula with a single explanatory
    variable like \code{~ z}. The observations in the model
    are ordered by the size of \code{z}. If set to \code{NULL} (the
    default) the observations are assumed to be ordered (e.g., a
    time series).}
 \item{prewhite}{logical or integer. Should the estimating functions
    be prewhitened? If \code{TRUE} or greater than 0 a VAR model of
    order \code{as.integer(prewhite)} is fitted via \code{ar} with
    method \code{"ols"} and \code{demean = FALSE}. The default is to
    use VAR(1) prewhitening.}
 \item{bw}{numeric. The bandwidth of the kernel (corresponds to the
    truncation lag). If set to \code{NULL} (the default) it is adaptively
    chosen by the function \code{bwAndrews}.}
 \item{kernel}{a character specifying the kernel used. All kernels used
    are described in Andrews (1991).}
 \item{approx}{a character specifying the approximation method if the
    bandwidth \code{bw} has to be chosen by \code{bwAndrews}.}    
 \item{diagnostics}{logical. Should additional model diagnostics be returned?
   See \code{\link{vcovHAC}} for details.}
 \item{sandwich}{logical. Should the sandwich estimator be computed?
    If set to \code{FALSE} only the middle matrix is returned.}
 \item{ar.method}{character. The \code{method} argument passed to
   \code{\link{ar}} for prewhitening (only, not for bandwidth selection).}
 \item{tol}{numeric. Weights that exceed \code{tol} are used for computing
   the covariance matrix, all other weights are treated as 0.}
 \item{data}{an optional data frame containing the variables in the \code{order.by} 
    model. By default the variables are taken from the environment which
    the function is called from.}
 \item{\dots}{further arguments passed to \code{bwAndrews}.}
 \item{weights}{numeric. A vector of weights used for weighting the estimated
   coefficients of the approximation model (as specified by \code{approx}). By
   default all weights are 1 except that for the intercept term (if there is more than
   one variable).}
}

\details{\code{kernHAC} is a convenience interface to \code{\link{vcovHAC}} using 
\code{weightsAndrews}: first a weights function is defined and then \code{vcovHAC}
is called.

The kernel weights underlying \code{weightsAndrews}
are directly accessible via the function \code{\link{kweights}} and require
the specification of the bandwidth parameter \code{bw}. If this is not specified
it can be chosen adaptively by the function \code{bwAndrews} (except for the
\code{"Truncated"} kernel). The automatic bandwidth selection is based on
an approximation of the estimating functions by either AR(1) or ARMA(1,1) processes.
To aggregate the estimated parameters from these approximations a weighted sum
is used. The \code{weights} in this aggregation are by default all equal to 1
except that corresponding to the intercept term which is set to 0 (unless there
is no other variable in the model) making the covariance matrix scale invariant.

Further details can be found in Andrews (1991).

The estimator of Newey & West (1987) can be obtained using the \code{"Bartlett"}
kernel.
}

\value{
\code{kernHAC} returns the same type of object as \code{\link{vcovHAC}}
which is typically just the covariance matrix.

\code{weightsAndrews} returns a vector of weights.

\code{bwAndrews} returns the selected bandwidth parameter.
}


\references{
  Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.

  Newey WK & West KD (1987),
  A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix.
  \emph{Econometrica}, \bold{55},
  703--708.
}

\seealso{\code{\link{vcovHAC}}, \code{\link{weightsLumley}},
         \code{\link{weave}}}

\examples{
curve(kweights(x, kernel = "Quadratic", normalize = TRUE),
      from = 0, to = 3.2, xlab = "x", ylab = "k(x)")
curve(kweights(x, kernel = "Bartlett", normalize = TRUE),
      from = 0, to = 3.2, col = 2, add = TRUE)
curve(kweights(x, kernel = "Parzen", normalize = TRUE),
      from = 0, to = 3.2, col = 3, add = TRUE)
curve(kweights(x, kernel = "Tukey", normalize = TRUE),
      from = 0, to = 3.2, col = 4, add = TRUE)
curve(kweights(x, kernel = "Truncated", normalize = TRUE),
      from = 0, to = 3.2, col = 5, add = TRUE)

x <- sin(1:100)
y <- 1 + x + rnorm(100)
fm <- lm(y ~ x)
kernHAC(fm)
vcov(fm)
}

\keyword{regression}
\keyword{ts}
