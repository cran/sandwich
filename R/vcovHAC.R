vcovHAC <- function(x, order.by = NULL, prewhite = FALSE,
  weights = weightsAndrews, diagnostics = FALSE, sandwich = TRUE,
  ar.method = "ols", data = list())
{
  prewhite <- as.integer(prewhite)

  umat <- estfun(x)[, , drop = FALSE]
  n.orig <- n <- nrow(umat)
  k <- ncol(umat)

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    index <- order(z)
  } else {
    index <- 1:n
  }
  umat <- umat[index, , drop = FALSE]

  if(prewhite > 0) {
    var.fit <- ar(umat, order.max = prewhite, demean = FALSE, aic = FALSE, method = ar.method)
    if(k > 1) D <- solve(diag(ncol(umat)) - apply(var.fit$ar, 2:3, sum))
      else D <- as.matrix(1/(1 - sum(var.fit$ar)))
    umat <- as.matrix(na.omit(var.fit$resid))
    n <- n - prewhite
  }

  if(is.function(weights))
    weights <- weights(x, order.by = order.by, prewhite = prewhite, ar.method = ar.method, data = data)

  if(length(weights) > n) {
    warning("more weights than observations, only first n used")
    weights <- weights[1:n]
  }
 
  utu <- 0.5 * crossprod(umat) * weights[1]
  wsum<-n*weights[1]/2
  w2sum<-n*weights[1]^2/2

  if(length(weights) > 1) {
    for (ii in 2:length(weights)) {
      utu <- utu + weights[ii] * crossprod(umat[1:(n-ii+1),,drop=FALSE], umat[ii:n,,drop=FALSE])
      wsum <- wsum + (n-ii+1) * weights[ii]
      w2sum <- w2sum + (n-ii+1) * weights[ii]^2
    }
  }

  utu <- utu + t(utu)
  ## Andrews: multiply with df n/(n-k)
  utu <- n.orig/(n.orig-k) * utu
  
  if(prewhite > 0) {
    utu <- crossprod(t(D), utu) %*% t(D)
  }
  
  
  wsum <- 2*wsum
  w2sum <- 2*w2sum
  bc <- n^2/(n^2 - wsum)
  df <- n^2/w2sum

  if(sandwich) {
    modelv <- summary(x)$cov.unscaled
    rval <- modelv %*% utu %*% modelv
  } else rval <- utu/n.orig

  if(diagnostics)  attr(rval, "diagnostics") <- list(bias.correction = bc, df = df)
  return(rval)
}

weightsAndrews <- function(x, order.by = NULL, bw = NULL,
  kernel = c("Quadratic Spectral", "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),
  prewhite = 1, ar.method = "ols", tol = 1e-7, data = list(), ...)
{
  kernel <- match.arg(kernel)
  if(is.null(bw))
    bw <- bwAndrews(x, order.by = order.by, kernel = kernel,
      prewhite = prewhite, data = data, ar.method = ar.method, ...)
  n <- length(residuals(x)) - as.integer(prewhite)
  
  weights <- kweights(0:(n-1)/bw, kernel = kernel)
  weights <- weights[1:max(which(abs(weights) > tol))]
  return(weights)
}

bwAndrews <- function(x, order.by = NULL, kernel = c("Quadratic Spectral", "Truncated",
  "Bartlett", "Parzen", "Tukey-Hanning"), approx = c("AR(1)", "ARMA(1,1)"),
  weights = NULL, prewhite = 1, ar.method = "ols", data = list())
{
  kernel <- match.arg(kernel)
  approx <- match.arg(approx)
  prewhite <- as.integer(prewhite)

  umat <- estfun(x)[,, drop = FALSE]
  n <- nrow(umat)
  k <- ncol(umat)

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    index <- order(z)
  } else {
    index <- 1:n
  }

  umat <- umat[index, , drop = FALSE]

  ## compute weights (try to set the intercept weight to 0)
  ## ignore at the moment: use 1 instead
  weights <- 1
  
  if(is.null(weights)) {
    weights <- rep(1, k)
    unames <- colnames(umat)
    if(!is.null(unames) && "(Intercept)" %in% unames)
      weights[which(unames == "(Intercept)")] <- 0
    else {
      res <- as.vector(residuals(x, "working"))
      weights[which(colSums((umat - res)^2) < 1e-16)] <- 0      
    }
  } else {
    weights <- rep(weights, length.out = k)
  }
  if(length(weights) < 2) weights <- 1

  if(prewhite > 0) {
    umat <- as.matrix(na.omit(ar(umat, order.max = prewhite,
      demean = FALSE, aic = FALSE, method = ar.method)$resid))
    n <- n - prewhite ##??
  }

  if(approx == "AR(1)") {
    fitAR1 <- function(x) {
      rval <-  ar(x, order.max = 1, aic = FALSE, method = "ols")
      rval <- c(rval$ar, sqrt(rval$var.pred))
      names(rval) <- c("rho", "sigma")
      return(rval)
    }

    ar.coef <- apply(umat, 2, fitAR1)

    denum <- sum(weights * (ar.coef["sigma",]/(1-ar.coef["rho",]))^4)
    alpha2 <- sum(weights * 4 * ar.coef["rho",]^2 * ar.coef["sigma",]^4/(1-ar.coef["rho",])^8) / denum
    alpha1 <- sum(weights * 4 * ar.coef["rho",]^2 * ar.coef["sigma",]^4/((1-ar.coef["rho",])^6 * (1+ar.coef["rho",])^2)) / denum

  } else {

    fitARMA11 <- function(x) {
      rval <-  arima(x, order = c(1, 0, 1), include.mean = FALSE)
      rval <- c(rval$coef, sqrt(rval$sigma2))
      names(rval) <- c("rho", "psi", "sigma")
      return(rval)
    }

    arma.coef <- apply(umat, 2, fitARMA11)

    denum <- sum(weights * ((1 + arma.coef["psi",]) * arma.coef["sigma",]/(1-arma.coef["rho",]))^4)
    alpha2 <- sum(weights * 4 * ((1 + arma.coef["rho",] * arma.coef["psi",]) * (
                                  arma.coef["rho",] + arma.coef["psi",]))^2 * arma.coef["sigma",]^4/
				 (1-arma.coef["rho",])^8) / denum
    alpha1 <- sum(weights * 4 * ((1 + arma.coef["rho",] * arma.coef["psi",]) * (
                                  arma.coef["rho",] + arma.coef["psi",]))^2 * arma.coef["sigma",]^4/
                                 ((1-arma.coef["rho",])^6 * (1+arma.coef["rho",])^2)) / denum
  }

  rval <- switch(kernel,
    "Truncated"          = {stop("no automatic bandwidth selection available for `Truncated' kernel, use weightsLumley() instead")},
    "Bartlett"           = {1.1447 * (n * alpha1)^(1/3)},
    "Parzen"             = {2.6614 * (n * alpha2)^(1/5)},
    "Tukey-Hanning"      = {1.7462 * (n * alpha2)^(1/5)},
    "Quadratic Spectral" = {1.3221 * (n * alpha2)^(1/5)})

  return(rval)  
}

weightsLumley <- function(x, order.by = NULL, C = NULL,
  method = c("truncate", "smooth"), acf = isoacf, data = list(), ...)
{
  method <- match.arg(method)
  res <- residuals(x, "response")
  n <- length(res)

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    index <- order(z)
  } else {
    index <- 1:n
  }
  res <- res[index]

  rhohat <- acf(res)

  switch(method,
  "truncate" = {
    if(is.null(C)) C <- 4
    lag <- max((1:length(rhohat))[rhohat^2*n > C])
    weights <- rep(1, lag)
  },
  "smooth" = {
    if(is.null(C)) C <- 1
    weights <- C * n * rhohat^2
    weights <- ifelse(weights > 1, 1, weights)
    weights <- weights[1:max(which(abs(weights) > 1e-7))]
  })
  
  return(weights)
}



kernHAC <- function(x, order.by = NULL, prewhite = 1, bw = NULL,
  kernel = c("Quadratic Spectral", "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),
  approx = c("AR(1)", "ARMA(1,1)"), diagnostics = FALSE, sandwich = TRUE,
  ar.method = "ols", tol = 1e-7, data = list(), ...)
{
  myweights <- function(x, order.by = NULL, prewhite = FALSE, ar.method = "ols", data = list())
    weightsAndrews(x, order.by = order.by, prewhite = prewhite, bw = bw,
                   kernel = kernel, approx = approx, ar.method = ar.method, tol = 1e-7, data = data, ...)
  vcovHAC(x, order.by = order.by, prewhite = prewhite,
    weights = myweights, diagnostics = diagnostics, sandwich = sandwich,
    ar.method = ar.method, data = data)
}

weave <- function(x, order.by = NULL, prewhite = FALSE, C = NULL,
  method = c("truncate", "smooth"), acf = isoacf, diagnostics = FALSE, 
  sandwich = TRUE, data = list(), ...)
{
  myweights <- function(x, order.by = NULL, prewhite = FALSE, data = list(), ...)
    weightsLumley(x, order.by = order.by, prewhite = prewhite, C = C,
                   method = method, acf = acf, data = data)
  vcovHAC(x, order.by = order.by, prewhite = prewhite,
    weights = myweights, diagnostics = diagnostics, sandwich = sandwich,
    data = data)
}
