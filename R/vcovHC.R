vcovHC <- function(x, 
  type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4"),
  omega = NULL, ...)
{
  if(is.matrix(x$x))
    X <- x$x
  else {
    mf <- model.frame(x)
    X <- model.matrix(terms(x), mf)    
  }
  res <- residuals(x)
  n <- nrow(X)
  k <- ncol(X)
  
  x.sum <- summary(x)  
  Q1 <- x.sum$cov.unscaled
  sigma2 <- x.sum$sigma^2
  diaghat <- hat(X)
  type <- match.arg(type)
  if(type == "HC") type <- "HC0"

  V <- NULL
  if(is.null(omega)) {
    switch(type,
      "const" = { omega <- function(residuals, diaghat, df) rep(1, length(residuals)) * sum(residuals^2)/df
                  V <- sigma2 * Q1 },
      "HC0" = { omega <- function(residuals, diaghat, df) residuals^2 },
      "HC1" = { omega <- function(residuals, diaghat, df) residuals^2 * length(residuals)/df },
      "HC2" = { omega <- function(residuals, diaghat, df) residuals^2 / (1 - diaghat) },
      "HC3" = { omega <- function(residuals, diaghat, df) residuals^2 / (1 - diaghat)^2 },
      "HC4" = { omega <- function(residuals, diaghat, df) residuals^2 / (1 - diaghat)^pmin(4, length(residuals) * diaghat/as.integer(round(sum(diaghat), digits = 0))) })
  }
  if(is.null(V)) {
    if(is.function(omega)) omega <- omega(res, diaghat, x$df.residual)
    VX <- sqrt(omega) * X
    V <- crossprod(crossprod(t(VX), Q1))
  }
  return(V)
}
