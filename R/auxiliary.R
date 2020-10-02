nobs0	<- function(x, ...) {
  nobs1 <- if("stats4" %in% loadedNamespaces()) stats4::nobs else stats::nobs
  nobs2 <- function(x, ...) NROW(residuals(x, ...))
  rval <- try(nobs1(x, ...), silent = TRUE)
  if(inherits(rval, "try-error") | is.null(rval)) rval <- nobs2(x, ...)
  return(rval)
}

vcov0	<- if(!is.null(vcov)) vcov else {
  if("stats4" %in% loadedNamespaces()) stats4::vcov else stats::vcov
}

kweights <- function(x, kernel = c("Truncated", "Bartlett", "Parzen",
                     "Tukey-Hanning", "Quadratic Spectral"), normalize = FALSE)
{
  kernel <- match.arg(kernel,
    c("Truncated", "Bartlett", "Parzen", "Tukey-Hanning", "Quadratic Spectral"))
  if(normalize) {
    ca <- switch(kernel,  
      "Truncated" = 2,
      "Bartlett" = 2/3,
      "Parzen" = .539285,
      "Tukey-Hanning" = 3/4,
      "Quadratic Spectral" = 1)
  } else {
    ca <- 1
  }

  switch(kernel,  
    "Truncated" = {
      ifelse(ca * abs(x) > 1, 0, 1)
    },
    "Bartlett" = {
      ifelse(ca * abs(x) > 1, 0, 1 - abs(ca * x))
    },
    "Parzen" = { 
      ifelse(ca * abs(x) > 1, 0, ifelse(ca * abs(x) < 0.5,
        1 - 6 * (ca * x)^2 + 6 * abs(ca * x)^3,
	2 * (1 - abs(ca * x))^3))
    },
    "Tukey-Hanning" = {
      ifelse(ca * abs(x) > 1, 0, (1 + cos(pi * ca * x))/2)
    },
    "Quadratic Spectral" = {
      qs <- function(x) {
        y <- 6 * pi * x/5
	3 * (1/y)^2 * (sin(y)/y - cos(y))
      }
      w <- qs(x)
      if(length(ix <- which(abs(x) < 1e-3)) > 0L) {
        cf <- 1e6 * log(qs(1e-3))
        w[ix] <- exp(cf * x[ix]^2)
      }
      w
    }
  )
}

isoacf <- function(x, lagmax = NULL, weave1 = FALSE)
{
  acfWeave <- function(x, lag = trunc(5*sqrt(length(x))))
  {
    x <- x - mean(x)
    # n <- length(x)
    autocov <- function(ii, xx)
      cov(xx[1:(length(xx)-ii+1)],xx[ii:length(xx)])
    covs <- sapply(2:lag, autocov, xx = x)
    covs/var(x)
  }

  if(weave1) {
    lagmax <- trunc(5*sqrt(length(x)))
    lagmax <- min(length(x) - 1, lagmax)
    covs <- acfWeave(x, lag = lagmax)
    isocovs <- pava.blocks(c(covs,0), c((length(x)-1):(length(x)-length(covs)),
      .Machine$double.xmax), up = FALSE)
    rval <- c(1, rep(isocovs$x, isocovs$blocks))    
  } else {
    lagmax <- length(x) - 1
    lagmax <- min(length(x) - 1, lagmax)
    covs <- as.vector(acf(x, lag.max = lagmax -1, plot = FALSE)$acf)[-1]
    rval <- c(1, -isoreg(1:(length(covs)+1), c(-covs, 0))$yf)
  }
  return(rval)
}

pava.blocks <- function(x, w = rep(1, length(x)),
  b = rep(1, length(x)), up = TRUE)
{
  lasti <- 1
  if(length(x) == 1) rval <- list(x = x, blocks = b,increasing = up)
  else {
    for(i in 2:length(x)) {
      if(x[i] <= x[lasti] & up){
        wtotal <- w[lasti]+w[i]
        x[lasti] <- (x[lasti]*w[lasti]+x[i]*w[i])/wtotal
        w[lasti] <- wtotal
        b[lasti] <- b[i]+b[lasti]
        b[i] <- 0
      } else if(x[i] <= x[lasti] & !up) {
        lasti <- i
      } else if(x[i] > x[lasti] & !up) {
        wtotal <- w[lasti]+w[i]
        x[lasti] <- (x[lasti]*w[lasti]+x[i]*w[i])/wtotal
        w[lasti] <- wtotal
        b[lasti] <- b[i]+b[lasti]
        b[i] <- 0
      } else if(x[i] > x[lasti] & up) {
        lasti <- i
      }
    }
  
    if(any(b == 0)) {
      subset <- (b > 0)
      rval <- pava.blocks(x[subset],w[subset],b[subset],up)
    } else
      rval <- list(x = x,blocks = b,increasing = up)
  }
  return(rval)  
}
