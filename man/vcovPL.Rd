\name{vcovPL}
\alias{vcovPL}
\alias{meatPL}
\encoding{UTF-8}

\title{Clustered Covariance Matrix Estimation for Panel Data}

\description{
  Estimation of sandwich covariances a la Newey-West (1987)
  and Driscoll and Kraay (1998) for panel data.
}

\usage{
vcovPL(x, cluster = NULL, order.by = NULL,
  kernel = "Bartlett", sandwich = TRUE, fix = FALSE, \dots)

meatPL(x, cluster = NULL, order.by = NULL,
  kernel = "Bartlett", lag = "NW1987", bw = NULL,
  adjust = TRUE, aggregate = TRUE, \dots)
}

\arguments{
  \item{x}{a fitted model object.}
  \item{cluster}{a single variable indicating the clustering of observations,
    or a \code{list} (or \code{data.frame}) of one or two variables, or a
    formula specifying which one ore two variables from the fitted model should
    be used (see examples). In case two variables are specified, the second variable
    is assumed to provide the time ordering (instead of using the argument
    \code{order.by}).
    By default (\code{cluster = NULL}), either \code{attr(x, "cluster")} is used
    (if any) or otherwise every observation is assumed to be its own cluster.}
  \item{order.by}{a variable, list/data.frame, or formula indicating the
    aggregation within time periods. By default \code{attr(x, "order.by")} is used
    (if any) or specified through the second variable in \code{cluster} (see above).
    If neither is available,  observations within clusters are assumed to be ordered.}
  \item{kernel}{a character specifying the kernel used. All kernels
    described in Andrews (1991) are supported, see \code{\link{kweights}}.}
  \item{lag}{character or numeric, indicating the lag length used.
    Three rules of thumb (\code{"max"} or equivalently \code{"P2009"},
    \code{"NW1987"}, or \code{"NW1994"}) can be specified, or a numeric
    number of lags can be specified directly. By default, \code{"NW1987"} is used.}
  \item{bw}{numeric. The bandwidth of the kernel which by default corresponds
    to \code{lag + 1}. Only one of \code{lag} and \code{bw} should be
    used.}
  \item{sandwich}{logical. Should the sandwich estimator be computed?
    If set to \code{FALSE} only the meat matrix is returned.}
  \item{fix}{logical. Should the covariance matrix be fixed to be
    positive semi-definite in case it is not?}
  \item{adjust}{logical. Should a finite sample adjustment be made? This
    amounts to multiplication with \eqn{n/(n - k)} where \eqn{n} is the
    number of observations and \eqn{k} is the number of estimated parameters.}
  \item{aggregate}{logical. Should the \code{estfun} be aggregated within each
    time period (yielding Driscoll and Kraay 1998) or not (restricting cross-sectional
    and cross-serial correlation to zero, yielding panel Newey-West)?}
  \item{\dots}{arguments passed to the \code{metaPL} or \code{estfun} function,
    respectively.}
}

\details{
  \code{vcovPL} is a function for estimating the Newey-West (1987) and
  Driscoll and Kraay (1998) covariance matrix.
  Driscoll and Kraay (1998) apply a Newey-West type correction to the
  sequence of cross-sectional averages of the moment conditions (see
  Hoechle (2007)). For large \eqn{T} (and regardless of the length of the
  cross-sectional dimension), the Driscoll and Kraay (1998)
  standard errors are robust to general forms of cross-sectional and
  serial correlation (Hoechle (2007)).
  The Newey-West (1987) covariance matrix restricts the Driscoll and
  Kraay (1998) covariance matrix to no cross-sectional correlation.
  
  The function \code{meatPL} is the work horse for estimating
  the meat of Newey-West (1987) and Driscoll and Kraay (1998)
  covariance matrix estimators. \code{vcovPL} is a wrapper calling
  \code{\link{sandwich}} and \code{\link{bread}} (Zeileis 2006).

  Default lag length is the \code{"NW1987"}. 
  For \code{lag = "NW1987"}, the lag length is chosen from the heuristic
  \eqn{floor[T^{(1/4)}]}. More details on lag length selection in Hoechle (2007). 
  For \code{lag = "NW1994"}, the lag length is taken from the first step
  of Newey and West's (1994) plug-in procedure.
  
  The \code{cluster}/\code{order.by} specification can be made in a number of ways:
  Either both can be a single variable or \code{cluster} can be a 
  \code{list}/\code{data.frame} of two variables.
  If \code{\link[stats]{expand.model.frame}} works for the model object \code{x},
  the \code{cluster} (and potentially additionally \code{order.by}) can also be
  a \code{formula}. By default (\code{cluster = NULL, order.by = NULL}),
  \code{attr(x, "cluster")} and \code{attr(x, "order.by")} are checked and
  used if available. If not, every observation is assumed to be its own cluster,
  and observations within clusters are assumed to be ordered accordingly.
  If the number of observations in the model \code{x} is smaller than in the
  original \code{data} due to \code{NA} processing, then the same \code{NA} processing
  can be applied to \code{cluster} if necessary (and \code{x$na.action} being
  available).
}

\value{
  A matrix containing the covariance matrix estimate.
}

\references{
Andrews DWK (1991).
  \dQuote{Heteroscedasticity and Autocorrelation Consistent Covariance Matrix Estimation},
  \emph{Econometrica}, 817--858.

Driscoll JC &  Kraay AC (1998).
  \dQuote{Consistent Covariance Matrix Estimation with Spatially Dependent Panel Data},
  \emph{The Review of Economics and Statistics}, \bold{80}(4), 549--560.

Hoechle D (2007).
  \dQuote{Robust Standard Errors for Panel Regressions with Cross-Sectional Dependence},
  \emph{Stata Journal}, \bold{7}(3), 281--312.

Newey WK & West KD (1987).
  \dQuote{Hypothesis Testing with Efficient Method of Moments Estimation},
  \emph{International Economic Review}, 777-787.

Newey WK & West KD (1994).
  \dQuote{Automatic Lag Selection in Covariance Matrix Estimation},
  \emph{The Review of Economic Studies}, \bold{61}(4), 631--653.
  
White H (1980).
  \dQuote{A Heteroskedasticity-Consistent Covariance Matrix Estimator and a Direct Test for Heteroskedasticity},
  \emph{Econometrica}, 817--838.
  \doi{10.2307/1912934}
  
Zeileis A (2004).
  \dQuote{Econometric Computing with HC and HAC Covariance Matrix Estimator},
  \emph{Journal of Statistical Software}, \bold{11}(10), 1--17.
  \doi{10.18637/jss.v011.i10}

Zeileis A (2006).
  \dQuote{Object-Oriented Computation of Sandwich Estimators},
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  \doi{10.18637/jss.v016.i09}

Zeileis A, Köll S, Graham N (2020).
  \dQuote{Various Versatile Variances: An Object-Oriented Implementation of Clustered Covariances in R.}
  \emph{Journal of Statistical Software}, \bold{95}(1), 1--36.
  \doi{10.18637/jss.v095.i01}
}

\seealso{\code{\link{vcovCL}}, \code{\link{vcovPC}}}

\examples{
## Petersen's data
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

## Driscoll and Kraay standard errors
## lag length set to: T - 1 (maximum lag length)
## as proposed by Petersen (2009)
sqrt(diag(vcovPL(m, cluster = ~ firm + year, lag = "max", adjust = FALSE)))

## lag length set to: floor(4 * (T / 100)^(2/9))
## rule of thumb proposed by Hoechle (2007) based on Newey & West (1994)
sqrt(diag(vcovPL(m, cluster = ~ firm + year, lag = "NW1994")))

## lag length set to: floor(T^(1/4))
## rule of thumb based on Newey & West (1987)
sqrt(diag(vcovPL(m, cluster = ~ firm + year, lag = "NW1987")))

## the following specifications of cluster/order.by are equivalent
vcovPL(m, cluster = ~ firm + year)
vcovPL(m, cluster = PetersenCL[, c("firm", "year")])
vcovPL(m, cluster = ~ firm, order.by = ~ year)
vcovPL(m, cluster = PetersenCL$firm, order.by = PetersenCL$year)

## these are also the same when observations within each
## cluster are already ordered
vcovPL(m, cluster = ~ firm)
vcovPL(m, cluster = PetersenCL$firm)
}

\keyword{regression}
