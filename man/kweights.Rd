\name{kweights}
\alias{kweights}
\title{Kernel Weights}
\description{
 Kernel weights for kernel-based heteroskedasticity
 and autocorrelation consistent (HAC) covariance matrix estimators
 as introduced by Andrews (1991).
}

\usage{
kweights(x, kernel = c("Truncated", "Bartlett", "Parzen",
  "Tukey-Hanning", "Quadratic Spectral"), normalize = FALSE)
}

\arguments{
 \item{x}{numeric.}
 \item{kernel}{a character specifying the kernel used. All kernels used
    are described in Andrews (1991).}
 \item{normalize}{logical. If set to \code{TRUE} the kernels are
   normalized as described in Andrews (1991).}    
}

\value{
Value of the kernel function at \code{x}.
}

\references{
  Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.
}

\seealso{\code{\link{kernHAC}}, \code{\link{weightsAndrews}}}

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
}

\keyword{regression}
\keyword{ts}
