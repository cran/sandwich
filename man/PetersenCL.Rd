\name{PetersenCL}
\alias{PetersenCL}

\title{Petersen's Simulated Data for Assessing Clustered Standard Errors}

\description{
  Artificial balanced panel data set from Petersen (2009) for
  illustrating and benchmarking clustered standard errors.
}

\usage{data("PetersenCL")}

\format{
A data frame containing 5000 observations on 4 variables.
\describe{
  \item{firm}{integer. Firm identifier (500 firms).}
  \item{year}{integer. Time variable (10 years per firm).}
  \item{x}{numeric. Independent regressor variable.}
  \item{y}{numeric. Dependent response variable.}
}
}

\details{
This simulated data set was created to illustrate and benchmark clustered standard errors.
The residual and the regressor variable both contain a firm effect, but no year effect.
Thus, standard errors clustered by firm are different from the OLS standard errors
and similarly double-clustered standard errors (by firm and year) are different from
the standard errors clustered by year.
}

\source{
\url{https://www.kellogg.northwestern.edu/faculty/petersen/htm/papers/se/test_data.htm}
}

\references{
Petersen MA (2009).
  \dQuote{Estimating Standard Errors in Finance Panel Data Sets: Comparing Approaches},
  \emph{The Review of Financial Studies}, \bold{22}(1), 435--480.
  \doi{10.1093/rfs/hhn053}
}

\keyword{datasets}
