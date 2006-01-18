estfun <- function(x, ...)
{
  UseMethod("estfun")
}

estfun.lm <- function(x, ...)
{
  if (is.matrix(x$x))
    xmat <- x$x
  else {
    mf <- model.frame(x)
    xmat <- model.matrix(terms(x), mf)    
  }
  wts <- if(!is.null(x$weights)) x$weights else 1
  res <- residuals(x)
  rval <- as.vector(res) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  if(is.zoo(res)) rval <- zoo(rval, index(res), attr(res, "frequency"))
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  return(rval)
}

estfun.glm <- function(x, ...)
{
  if (is.matrix(x$x))
    xmat <- x$x
  else {
    mf <- model.frame(x)
    xmat <- model.matrix(terms(x), mf)    
  }
  rval <- as.vector(residuals(x, "working")) * x$weights * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  res <- residuals(x, type = "pearson")
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  if(is.zoo(res)) rval <- zoo(rval, index(res), attr(res, "frequency"))
  return(rval)
}

estfun.rlm <- function(x, ...)
{
  if (is.matrix(x$x)) 
      xmat <- x$x
  else {
      mf <- model.frame(x)
      xmat <- model.matrix(terms(x), mf)
  }
  if (!is.null(x$weights)) 
      wts <- x$weights
  else wts <- 1
  res <- residuals(x)
  psi <- function(z) x$psi(z) * z
  rval <- as.vector(psi(res/x$s)) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  if(is.zoo(res)) rval <- zoo(rval, index(res), attr(res, "frequency"))
  return(rval)
}

## check: coxph and survreg

estfun.coxph <- function(x, ...)
{
  stopifnot(require("survival"))
  residuals(x, type = "score", ...)
}

estfun.survreg <- function(x, ...)
{
  stopifnot(require("survival"))
  if (is.matrix(x$x))
    xmat <- x$x
  else {
    mf <- model.frame(x)
    xmat <- model.matrix(terms(x), mf)    
  }
  wts <- if(!is.null(x$weights)) x$weights else 1
  res <- residuals(x, type = "matrix")
  rval <- as.vector(res[,"dg"]) * wts * xmat
  if(NROW(x$var) > length(coef(x))) {
    rval <- cbind(rval, res[,"ds"])
    colnames(rval)[NCOL(rval)] <- "Log(scale)"
  }
  attr(rval, "assign") <- NULL
  
  return(rval)
}

estfun.nls <- function(x, ...)
{
  rval <- as.vector(x$m$resid()) * x$m$gradient()
  colnames(rval) <- names(coef(x))
  rval
}
