vcovHC <- function(x, order.by = NULL, data = list(),
  type = c("HC2", "const", "HC", "HC1", "HC3"))
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

  if(!is.null(order.by))
  {
    z <- model.matrix(order.by, data = data)
    z <- as.vector(z[,ncol(z)])
    index <- order(z)
  } else {
    index <- 1:n
  }
  X <- X[index, , drop = FALSE]
  res <- res[index]
    
  Q1 <- summary(x)$cov.unscaled
  sigma2 <- var(res)*(n-1)/(n-k)
  type <- match.arg(type)

  if( type == "const") {
    V <- sigma2 * Q1
  } else {
    if(type == "HC2")
    {
      diaghat <- 1 - diag(X %*% Q1 %*% t(X))
      res <- res/sqrt(diaghat)
    }
    if(type == "HC3")
    {
      diaghat <- 1 - diag(X %*% Q1 %*% t(X))
      res <- res/diaghat
      Xu <- as.vector(t(X) %*% res)
    }
    VX <- res * X
    if(type %in% c("HC", "HC1", "HC2")) { V <- crossprod(crossprod(t(VX), Q1)) }
    if(type == "HC1") {V <- V * (n/(n-k))}
    if(type == "HC3") {V <- Q1 %*% (crossprod(VX) - (outer(Xu,Xu) /n)) %*% Q1 * (n-1)/n}
  }
  return(V)
}

