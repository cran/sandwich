sandwich <- function(x, bread. = bread, meat. = meat, ...)
{
  if(is.function(bread.)) bread. <- bread.(x)
  if(is.function(meat.)) meat. <- meat.(x, ...)
  n <- NROW(estfun(x))
  return(1/n * (bread. %*% meat. %*% bread.))
}

bread <- function(x, ...)
{
  UseMethod("bread")
}

bread.lm <- bread.nls <- function(x, ...)
{
  sx <- summary(x)
  return(sx$cov.unscaled * as.vector(sum(sx$df[1:2])))
}

bread.survreg <- function(x, ...)
  length(x$linear.predictors) * x$var

bread.gam <- function(x, ...)
{
  sx <- summary(x)
  sx$cov.unscaled * sx$n
}
  
bread.coxph <- function(x, ...)
{
  rval <- x$var * x$n
  dimnames(rval) <- list(names(coef(x)), names(coef(x)))
  return(rval)

}

meat <- function(x, adjust = FALSE, ...)
{
  psi <- estfun(x, ...)
  k <- NCOL(psi)
  n <- NROW(psi)
  rval <- crossprod(as.matrix(psi))/n
  if(adjust) rval <- n/(n-k) * rval
  rownames(rval) <- colnames(rval) <- colnames(psi)
  return(rval)
}
