bread <- function(x, ...)
{
  UseMethod("bread")
}

bread.lm <- function(x, ...)
{
  sx <- summary.lm(x)
  return(sx$cov.unscaled * as.vector(sum(sx$df[1:2])))
}

bread.glm <- function(x, ...)
{
  sx <- summary(x)
  wres <- as.vector(residuals(x, "working")) * weights(x, "working")
  dispersion <- if(substr(x$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
    else sum(wres^2)/sum(weights(x, "working"))
  return(sx$cov.unscaled * as.vector(sum(sx$df[1:2])) * dispersion)
}

bread.nls <- function(x, ...)
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

bread.hurdle <- function(x, ...) {
  x$vcov * x$n
}

bread.zeroinfl <- function(x, ...) {
  x$vcov * x$n
}
