\name{vcovPC}
\alias{vcovPC}
\alias{meatPC}
\encoding{UTF-8}

\title{Panel-Corrected Covariance Matrix Estimation}

\description{
  Estimation of sandwich covariances a la Beck and Katz (1995) for panel data.
}

\usage{
vcovPC(x, cluster = NULL, order.by = NULL,
  pairwise = FALSE, sandwich = TRUE, fix = FALSE, \dots)

meatPC(x, cluster = NULL, order.by = NULL,
  pairwise = FALSE, kronecker = TRUE, \dots)
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
  \item{pairwise}{logical. For unbalanced panels. Indicating whether
  the meat should be estimated pair- or casewise.}
  \item{sandwich}{logical. Should the sandwich estimator be computed?
    If set to \code{FALSE} only the meat matrix is returned.}
  \item{fix}{logical. Should the covariance matrix be fixed to be
    positive semi-definite in case it is not?}  
  \item{kronecker}{logical. Calculate the meat via the
  Kronecker-product, shortening the computation time for small
  matrices. For large matrices, set \code{kronecker = FALSE}.}
  \item{\dots}{arguments passed to the \code{meatPC} or \code{estfun} function,
    respectively.}
}

\details{
  \code{vcovPC} is a function for estimating Beck and Katz (1995)
  panel-corrected covariance matrix.
 
  The function \code{meatPC} is the work horse for estimating
  the meat of Beck and Katz (1995) covariance matrix estimators.
  \code{vcovPC} is a wrapper calling
  \code{\link[sandwich]{sandwich}} and \code{\link[sandwich]{bread}}
  (Zeileis 2006).

  Following Bailey and Katz (2011), there are two alternatives to
  estimate the meat for unbalanced panels.
  For \code{pairwise = FALSE}, a balanced subset of the panel is used,
  whereas for \code{pairwise = TRUE}, a pairwise balanced sample is
  employed.

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
Bailey D, Katz JN (2011).
  \dQuote{Implementing Panel-Corrected Standard Errors in R: The pcse Package},
  \emph{Journal of Statistical Software, Code Snippets}, \bold{42}(1), 1--11.
  \doi{10.18637/jss.v042.c01}
  
Beck N, Katz JN (1995).
  \dQuote{What To Do (and Not To Do) with Time-Series-Cross-Section Data in Comparative Politics},
  \emph{American Political Science Review}, \bold{89}(3), 634--647.
  \doi{10.2307/2082979}
  
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

\seealso{\code{\link{vcovCL}}, \code{\link{vcovPL}}}

\examples{
## Petersen's data
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

## Beck and Katz (1995) standard errors
## balanced panel
sqrt(diag(vcovPC(m, cluster = ~ firm + year)))

## unbalanced panel
PU <- subset(PetersenCL, !(firm == 1 & year == 10))
pu_lm <- lm(y ~ x, data = PU)
sqrt(diag(vcovPC(pu_lm, cluster = ~ firm + year, pairwise = TRUE)))
sqrt(diag(vcovPC(pu_lm, cluster = ~ firm + year, pairwise = FALSE)))

\donttest{
## the following specifications of cluster/order.by are equivalent
vcovPC(m, cluster = ~ firm + year)
vcovPC(m, cluster = PetersenCL[, c("firm", "year")])
vcovPC(m, cluster = ~ firm, order.by = ~ year)
vcovPC(m, cluster = PetersenCL$firm, order.by = PetersenCL$year)

## these are also the same when observations within each
## cluster are already ordered
vcovPC(m, cluster = ~ firm)
vcovPC(m, cluster = PetersenCL$firm)
}
}

\keyword{regression}
