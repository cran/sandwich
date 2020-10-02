\name{vcovCL}
\alias{vcovCL}
\alias{meatCL}
\encoding{UTF-8}

\title{Clustered Covariance Matrix Estimation}

\description{
  Estimation of one-way and multi-way clustered
  covariance matrices using an object-oriented approach.
}

\usage{
vcovCL(x, cluster = NULL, type = NULL, sandwich = TRUE, fix = FALSE, \dots)
meatCL(x, cluster = NULL, type = NULL, cadjust = TRUE, multi0 = FALSE, \dots)
}

\arguments{
  \item{x}{a fitted model object.}
  \item{cluster}{a variable indicating the clustering of observations,
    a \code{list} (or \code{data.frame}) thereof, or a formula specifying
    which variables from the fitted model should be used (see examples).
    By default (\code{cluster = NULL}), either \code{attr(x, "cluster")} is used
    (if any) or otherwise every observation is assumed to be its own cluster.}
  \item{type}{a character string specifying the estimation type (HC0--HC3). 
    The default is to use \code{"HC1"} for \code{lm} objects and
    \code{"HC0"} otherwise.}
  \item{sandwich}{logical. Should the sandwich estimator be computed?
    If set to \code{FALSE} only the meat matrix is returned.}
  \item{fix}{logical. Should the covariance matrix be fixed to be
    positive semi-definite in case it is not?}
  \item{cadjust}{logical. Should a cluster adjustment be applied?}
  \item{multi0}{logical. Should the HC0 estimate be used for
    the final adjustment in multi-way clustered covariances?}
  \item{\dots}{arguments passed to \code{meatCL}.}
}

\details{
  Clustered sandwich estimators are used to adjust inference when errors
  are correlated within (but not between) clusters. \code{vcovCL} allows
  for clustering in arbitrary many cluster dimensions (e.g., firm, time, industry), given all
  dimensions have enough clusters (for more details, see Cameron et al. 2011).
  If each observation is its own cluster, the clustered sandwich
  collapses to the basic sandwich covariance.
  
  The function \code{meatCL} is the work horse for estimating
  the meat of clustered sandwich estimators. \code{vcovCL} is a wrapper calling
  \code{\link{sandwich}} and \code{\link{bread}} (Zeileis 2006).
  \code{vcovCL} is applicable beyond \code{lm} or \code{glm} class objects.

  \code{\link{bread}} and \code{\link{meat}} matrices are multiplied to
  construct clustered sandwich estimators.
  The meat of a clustered sandwich estimator is the cross product of
  the clusterwise summed estimating functions. Instead of summing over
  all individuals, first sum over cluster.
 
  A two-way clustered sandwich estimator \eqn{M} (e.g., for cluster dimensions
  "firm" and "industry" or "id" and "time") is a linear combination of
  one-way clustered sandwich estimators for both dimensions
  (\eqn{M_{id}, M_{time}}) minus the
  clustered sandwich estimator, with clusters formed out of the
  intersection of both dimensions (\eqn{M_{id \cap time}}):
  \deqn{M = M_{id} + M_{time} - M_{id \cap time}}
  Additionally, each of the three terms can be weighted by the corresponding
  cluster bias adjustment factor (see below and Equation 20 in Zeileis et al. 2020).
  Instead of subtracting \eqn{M_{id \cap time}} as the last
  subtracted matrix, Ma (2014) suggests to subtract the basic HC0
  covariance matrix when only a single observation is in each
  intersection of \eqn{id}  and \eqn{time}.
  Set \code{multi0 = TRUE} to subtract the basic HC0 covariance matrix as
  the last subtracted matrix in multi-way clustering. For details,
  see also Petersen (2009) and Thompson (2011).
  
  With the \code{type} argument, HC0 to HC3 types of
  bias adjustment can be employed, following the terminology used by
  MacKinnon and White (1985) for heteroscedasticity corrections. HC0 applies no small sample bias adjustment.
  HC1 applies a degrees of freedom-based correction, \eqn{(n-1)/(n-k)} where \eqn{n} is the
  number of observations and \eqn{k} is the number of explanatory or predictor variables in the model.
  HC1 is the most commonly used approach, and is the default, though it is less effective
  than HC2 and HC3 when the number of clusters is relatively small (Cameron et al. 2008).
  HC2 and HC3 types of bias adjustment are geared towards the linear
  model, but they are also applicable for GLMs (see Bell and McCaffrey
  2002, and Kauermann and Carroll 2001, for details).
  A precondition for HC2 and HC3 types of bias adjustment is the availability
  of a hat matrix (or a weighted version therof for GLMs) and hence
  these two types are currently only implemented for \code{\link{lm}}
  and \code{\link{glm}} objects.
  
  The \code{cadjust} argument allows to
  switch the cluster bias adjustment factor \eqn{G/(G-1)} on and
  off (where \eqn{G} is the number of clusters in a cluster dimension \eqn{g})
  See Cameron et al. (2008) and Cameron et al. (2011) for more details about
  small-sample modifications.

  The \code{cluster} specification can be made in a number of ways: The \code{cluster}
  can be a single variable or a \code{list}/\code{data.frame} of multiple
  clustering variables. If \code{\link[stats]{expand.model.frame}} works
  for the model object \code{x}, the \code{cluster} can also be a \code{formula}.
  By default (\code{cluster = NULL}), \code{attr(x, "cluster")} is checked and
  used if available. If not, every observation is assumed to be its own cluster.
  If the number of observations in the model \code{x} is smaller than in the
  original \code{data} due to \code{NA} processing, then the same \code{NA} processing
  can be applied to \code{cluster} if necessary (and \code{x$na.action} being
  available).

  Cameron et al. (2011) observe that sometimes the covariance matrix is
  not positive-semidefinite and recommend to employ the eigendecomposition of the estimated
  covariance matrix, setting any negative eigenvalue(s) to zero. This fix
  is applied, if necessary, when \code{fix = TRUE} is specified.
}

\value{
  A matrix containing the covariance matrix estimate.
}

\references{
Bell RM, McCaffrey DF (2002).
  \dQuote{Bias Reduction in Standard Errors for Linear Regression with Multi-Stage Samples},
  \emph{Survey Methodology}, \bold{28}(2), 169--181.

Cameron AC, Gelbach JB, Miller DL (2008).
  \dQuote{Bootstrap-Based Improvements for Inference with Clustered Errors},
  \emph{The Review of Economics and Statistics}, \bold{90}(3),
  414--427.
  \doi{10.3386/t0344}
  
Cameron AC, Gelbach JB, Miller DL (2011).
  \dQuote{Robust Inference with Multiway Clustering},
  \emph{Journal of Business & Ecomomic Statistics}, \bold{29}(2),
  238--249.
  \doi{10.1198/jbes.2010.07136}
  
Kauermann G, Carroll RJ (2001).
  \dQuote{A Note on the Efficiency of Sandwich Covariance Matrix
    Estimation},
  \emph{Journal of the American Statistical Association},
  \bold{96}(456), 1387--1396.
  \doi{10.1198/016214501753382309}

Ma MS (2014).
  \dQuote{Are We Really Doing What We Think We Are Doing? A Note on
    Finite-Sample Estimates of Two-Way Cluster-Robust Standard Errors},
  \emph{Mimeo, Availlable at SSRN:}
   URL \url{https://www.ssrn.com/abstract=2420421}.

MacKinnon, JG, White, H (1985). 
  \dQuote{Some Heteroskedasticity-Consistent Covariance Matrix Estimators with Improved Finite Sample Properties}
  \emph{Journal of Econometrics}, \bold{29}(3), 305--325.
  \doi{10.1016/0304-4076(85)90158-7}
  
Petersen MA (2009).
  \dQuote{Estimating Standard Errors in Finance Panel Data Sets: Comparing Approaches},
  \emph{The Review of Financial Studies}, \bold{22}(1), 435--480.
  \doi{10.1093/rfs/hhn053}

Thompson SB (2011).
  \dQuote{Simple Formulas for Standard Errors That Cluster by Both Firm
    and Time},
  \emph{Journal of Financial Economics}, \bold{99}(1), 1--10.
  \doi{10.1016/j.jfineco.2010.08.016}
  
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

\seealso{\code{\link{vcovHC}}}

\examples{
## Petersen's data
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

## clustered covariances
## one-way
vcovCL(m, cluster = ~ firm)
vcovCL(m, cluster = PetersenCL$firm) ## same
## one-way with HC2
vcovCL(m, cluster = ~ firm, type = "HC2")
## two-way
vcovCL(m, cluster = ~ firm + year)
vcovCL(m, cluster = PetersenCL[, c("firm", "year")]) ## same

## comparison with cross-section sandwiches
## HC0
all.equal(sandwich(m), vcovCL(m, type = "HC0", cadjust = FALSE))
## HC2
all.equal(vcovHC(m, type = "HC2"), vcovCL(m, type = "HC2"))
## HC3
all.equal(vcovHC(m, type = "HC3"), vcovCL(m, type = "HC3"))

## Innovation data
data("InstInnovation", package = "sandwich")

## replication of one-way clustered standard errors for model 3, Table I
## and model 1, Table II in Berger et al. (2017), see ?InstInnovation

## count regression formula
f1 <- cites ~ institutions + log(capital/employment) + log(sales) + industry + year

## model 3, Table I: Poisson model
## one-way clustered standard errors
tab_I_3_pois <- glm(f1, data = InstInnovation, family = poisson)
vcov_pois <- vcovCL(tab_I_3_pois, InstInnovation$company)
sqrt(diag(vcov_pois))[2:4]

## coefficient tables
if(require("lmtest")) {
coeftest(tab_I_3_pois, vcov = vcov_pois)[2:4, ]
}

\dontrun{
## model 1, Table II: negative binomial hurdle model
## (requires "pscl" or alternatively "countreg" from R-Forge)
library("pscl")
library("lmtest")
tab_II_3_hurdle <- hurdle(f1, data = InstInnovation, dist = "negbin")
#  dist = "negbin", zero.dist = "negbin", separate = FALSE)
vcov_hurdle <- vcovCL(tab_II_3_hurdle, InstInnovation$company)
sqrt(diag(vcov_hurdle))[c(2:4, 149:151)]
coeftest(tab_II_3_hurdle, vcov = vcov_hurdle)[c(2:4, 149:151), ]
}
}

\keyword{regression}
