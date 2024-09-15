vcovHC <- function(x, ...) {
  UseMethod("vcovHC")
}

vcovHC.default <- function(x, 
  type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5"),
  omega = NULL, sandwich = TRUE, ...)
{
  type <- match.arg(type)
  rval <- meatHC(x, type = type, omega = omega)
  if(sandwich) rval <- sandwich(x, meat. = rval, ...)
  return(rval)
}

meatHC <- function(x, 
  type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5"),
  omega = NULL, ...)
{
  ## ensure that NAs are omitted
  if(is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"

  ## extract X and residual degrees of freedom
  X <- model.matrix(x)
  if(any(alias <- is.na(coef(x)))) X <- X[, !alias, drop = FALSE]
  attr(X, "assign") <- NULL
  n <- NROW(X)
  df <- n - NCOL(X)

  ## hat values
  diaghat <- try(hatvalues(x), silent = TRUE)
  
  ## the following might work, but "intercept" is also claimed for "coxph"
  ## res <- if(attr(terms(x), "intercept") > 0) estfun(x)[,1] else rowMeans(estfun(x)/X, na.rm = TRUE)
  ## hence better rely on
  ef <- estfun(x, ...)
  res <- rowMeans(ef/X, na.rm = TRUE)

  ## handle rows with just zeros
  all0 <- apply(abs(ef) < .Machine$double.eps, 1L, all)
  res[all0] <- 0
  ## in case of lm/glm and type = "const" re-obtain the working residuals
  if(any(all0) && substr(type, 1L, 1L) == "c") {
    if(inherits(x, "glm")) {
      res <- as.vector(residuals(x, "working")) * weights(x, "working")
      if(!(substr(x$family$family, 1L, 17L) %in% c("poisson", "binomial", "Negative Binomial"))) {
        res <- res * sum(weights(x, "working"), na.rm = TRUE) / sum(res^2, na.rm = TRUE)
      }
    } else if(inherits(x, "lm")) {
      res <- as.vector(residuals(x))
      if(!is.null(weights(x))) res <- res * weights(x)
    }
  }

  ## if omega not specified, set up using type
  if(is.null(omega)) {
    type <- match.arg(type)
    if(type == "HC") type <- "HC0"
    switch(type,
      "const" = { omega <- function(residuals, diaghat, df) rep(1, length(residuals)) * sum(residuals^2)/df },
      "HC0"   = { omega <- function(residuals, diaghat, df) residuals^2 },
      "HC1"   = { omega <- function(residuals, diaghat, df) residuals^2 * length(residuals)/df },
      "HC2"   = { omega <- function(residuals, diaghat, df) residuals^2 / (1 - diaghat) },
      "HC3"   = { omega <- function(residuals, diaghat, df) residuals^2 / (1 - diaghat)^2 },
      "HC4"   = { omega <- function(residuals, diaghat, df) {
        n <- length(residuals)
	p <- as.integer(round(sum(diaghat),  digits = 0))
	delta <- pmin(4, n * diaghat/p)
        residuals^2 / (1 - diaghat)^delta
      }},
      "HC4m"  = { omega <- function(residuals, diaghat, df) {
        gamma <- c(1.0, 1.5) ## as recommended by Cribari-Neto & Da Silva
        n <- length(residuals)
	p <- as.integer(round(sum(diaghat),  digits = 0))
	delta <- pmin(gamma[1], n * diaghat/p) + pmin(gamma[2], n * diaghat/p)
        residuals^2 / (1 - diaghat)^delta
      }},
      "HC5"   = { omega <- function(residuals, diaghat, df) {
        k <- 0.7 ## as recommended by Cribari-Neto et al.
        n <- length(residuals)
	p <- as.integer(round(sum(diaghat),  digits = 0))
	delta <- pmin(n * diaghat/p, pmax(4, n * k * max(diaghat)/p))
        residuals^2 / sqrt((1 - diaghat)^delta)
      }}
    )
    
    ## check hat values
    if(type != "const") {
      if(type %in% c("HC0", "HC1")) {
        idx <- if(inherits(diaghat, "try-error")) abs(res)/sqrt(mean(res^2)) < .Machine$double.eps^0.75 else diaghat > 1 - sqrt(.Machine$double.eps)
        msg <- "%s covariances become (close to) singular if hat values are (close to) 1 as for observation(s) %s"
      } else {
        if(inherits(diaghat, "try-error")) stop(sprintf("hatvalues() could not be extracted successfully but are needed for %s", type))
        idx <- diaghat > 1 - sqrt(.Machine$double.eps)
        msg <- "%s covariances are numerically unstable for hat values close to 1 (and undefined if exactly 1) as for observation(s) %s"
      }
      if(any(idx)) {
        idx <- if(is.null(rownames(X))) as.character(which(idx)) else rownames(X)[idx]
        if(length(idx) > 10L) idx <- c(idx[1L:10L], "...")
        warning(sprintf(msg, type, paste(idx, collapse = ", ")))
      }
    }
  }
  
  ## process omega
  if(is.function(omega)) omega <- omega(res, diaghat, df)
  rval <- sqrt(omega) * X
  rval <- crossprod(rval)/n

  return(rval)
}

vcovHC.mlm <- function(x, 
  type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5"),
  omega = NULL, sandwich = TRUE, ...)
{
  ## compute meat "by hand" because meatHC() can not deal with
  ## residual "matrices"
  
  ## ensure that NAs are omitted
  if(is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"

  ## regressors and df
  X <- model.matrix(x)
  attr(X, "assign") <- NULL
  if(any(alias <- apply(coef(x), 1, function(z) all(is.na(z))))) X <- X[, !alias, drop = FALSE]
  n <- NROW(X)
  df <- n - NCOL(X)
  
  ## hat values and residuals
  diaghat <- hatvalues(x)
  res <- residuals(x)

  ## if omega not specified, set up using type
  if(is.null(omega)) {
    type <- match.arg(type)
    if(type == "HC") type <- "HC0"
    switch(type,
      "const" = {
        warning("not implemented for 'mlm' objects, using vcov() instead")
	return(vcov(x))
       },
      "HC0"   = { omega <- function(residuals, diaghat, df) residuals^2 },
      "HC1"   = { omega <- function(residuals, diaghat, df) residuals^2 * length(residuals)/df },
      "HC2"   = { omega <- function(residuals, diaghat, df) residuals^2 / (1 - diaghat) },
      "HC3"   = { omega <- function(residuals, diaghat, df) residuals^2 / (1 - diaghat)^2 },
      "HC4"   = { omega <- function(residuals, diaghat, df) {
        n <- length(residuals)
	p <- as.integer(round(sum(diaghat),  digits = 0))
	delta <- pmin(4, n * diaghat/p)
        residuals^2 / (1 - diaghat)^delta
      }},
      "HC4m"  = { omega <- function(residuals, diaghat, df) {
        gamma <- c(1.0, 1.5)
        n <- length(residuals)
	p <- as.integer(round(sum(diaghat),  digits = 0))
	delta <- pmin(gamma[1], n * diaghat/p) + pmin(gamma[2], n * diaghat/p)
        residuals^2 / (1 - diaghat)^delta
      }},
      "HC5"   = { omega <- function(residuals, diaghat, df) {
        k <- 0.7
        n <- length(residuals)
	p <- as.integer(round(sum(diaghat),  digits = 0))
	delta <- pmin(n * diaghat/p, pmax(4, n * k * max(diaghat)/p))
        residuals^2 / sqrt((1 - diaghat)^delta)
      }}
    )

    ## check hat values (if necessary)
    if(type %in% c("HC2", "HC3", "HC4", "HC4m", "HC5")) {
      if(inherits(diaghat, "try-error")) stop(sprintf("hatvalues() could not be extracted successfully but are needed for %s", type))
      id <- which(diaghat > 1 - sqrt(.Machine$double.eps))
      if(length(id) > 0L) {
        id <- if(is.null(rownames(X))) as.character(id) else rownames(X)[id]
        if(length(id) > 10L) id <- c(id[1L:10L], "...")
        warning(sprintf("%s cannot be computed due to observations with hatvalue = 1: %s", type, paste(id, collapse = ", ")))
      }
    }
  }
  
  ## process omega
  omega <- apply(res, 2, function(r) omega(r, diaghat, df))

  ## compute crossproduct X' Omega X
  rval <- lapply(1:ncol(omega), function(i) as.vector(sign(res[,i]) * sqrt(omega[,i])) * X)
  rval <- do.call("cbind", rval)
  rval <- crossprod(rval)/n
  colnames(rval) <- rownames(rval) <- colnames(vcov(x))

  ## if necessary compute full sandwich
  if(sandwich) rval <- sandwich(x, meat. = rval, ...)
  return(rval)
}
