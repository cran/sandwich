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
  if(!is.null(x$weights)) wts <- x$weights
    else wts <- 1
  res <- residuals(x)
  rval <- as.vector(res) * wts * xmat
  attr(rval, "assign") <- NULL
  if(is.zoo(res)) rval <- zoo(rval, time(res))
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
  res <- residuals(x, "pearson")
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  if(is.zoo(res)) rval <- zoo(rval, time(res))
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
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  if(is.zoo(res)) rval <- zoo(rval, time(res))
  return(rval)
}

# estfun.coxph and estfun.survreg would be nice.
# This one seems pretty close: note, however, that the
# intercept score is missing.
#
# estfun.coxph <- function(x, ...)
# {
#   stopifnot(require(survival))
#   res <- residuals(x)
#   rval <- residuals(x, type = "score", ...)
#   if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
#   if(is.zoo(res)) rval <- zoo(rval, time(res))
#   return(rval)
# }
