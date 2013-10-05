\name{isoacf}
\alias{isoacf}
\alias{pava.blocks}
\title{Isotonic Autocorrelation Function}
\description{
 Autocorrelation function (forced to be decreasing by isotonic regression).
}

\usage{
isoacf(x, lagmax = NULL, weave1 = FALSE)
}

\arguments{
 \item{x}{numeric vector.}
 \item{lagmax}{numeric. The maximal lag of the autocorrelations.}
 \item{weave1}{logical. If set to \code{TRUE} \code{isoacf}
   uses the \code{acf.R} and \code{pava.blocks} function from the
   original \code{weave} package, otherwise R's own \code{acf} and
   \code{isoreg} functions are used.}
}

\details{
\code{isoacf} computes the autocorrelation function (ACF)
of \code{x} enforcing the ACF to be decreasing by isotonic regression.
See also Robertson et al. (1988).
}

\value{
\code{isoacf} returns a numeric vector containing the ACF.
}

\references{
  Lumley T & Heagerty P (1999),
  Weighted Empirical Adaptive Variance Estimators for Correlated Data Regression.
  \emph{Journal of the Royal Statistical Society B}, \bold{61},
  459--477.

  Robertson T, Wright FT, Dykstra RL (1988),
  \emph{Order Restricted Statistical Inference}.
  New York. Wiley.
}

\seealso{\code{\link{weave}}, \code{\link{weightsLumley}}}

\examples{
x <- filter(rnorm(100), 0.9, "recursive")
isoacf(x)
acf(x, plot = FALSE)$acf
}

\keyword{regression}
\keyword{ts}
