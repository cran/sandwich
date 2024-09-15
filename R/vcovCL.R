vcovCL <- function(x, cluster = NULL, type = NULL, sandwich = TRUE, fix = FALSE, ...)
{
  ## compute meat of sandwich
  rval <- meatCL(x, cluster = cluster, type = type, ...)

  ## full sandwich
  if(sandwich) rval <- sandwich(x, meat. = rval)

  ## check (and fix) if sandwich is not positive semi-definite
  if(fix && any((eig <- eigen(rval, symmetric = TRUE))$values < 0)) {
    eig$values <- pmax(eig$values, 0)
    rval[] <- crossprod(sqrt(eig$values) * t(eig$vectors))
  }
  return(rval)
}

meatCL <- function(x, cluster = NULL, type = NULL, cadjust = TRUE, multi0 = FALSE, ...)
{
  ## extract estimating functions / aka scores
  if (is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"
  ef <- estfun(x, ...)
  k <- NCOL(ef)
  n <- NROW(ef)

  ## set up return value with correct dimension and names
  rval <- matrix(0, nrow = k, ncol = k,
    dimnames = list(colnames(ef), colnames(ef)))

  ## cluster can either be supplied explicitly or
  ## be an attribute of the model...FIXME: other specifications?
  if (is.null(cluster)) cluster <- attr(x, "cluster")

  ## resort to cross-section if no clusters are supplied
  if (is.null(cluster)) cluster <- 1L:n

  ## collect 'cluster' variables in a data frame
  if(inherits(cluster, "formula")) {
    cluster_tmp <- if("Formula" %in% loadedNamespaces()) { ## FIXME to suppress potential warnings due to | in Formula
      suppressWarnings(expand.model.frame(x, cluster, na.expand = FALSE))
    } else {
      expand.model.frame(x, cluster, na.expand = FALSE)
    }
    cluster <- model.frame(cluster, cluster_tmp, na.action = na.pass)
  } else {
    cluster <- as.data.frame(cluster)
  }
  
  ## handle omitted or excluded observations
  if((n != NROW(cluster)) && !is.null(x$na.action) && (class(x$na.action) %in% c("exclude", "omit"))) {
    cluster <- cluster[-x$na.action, , drop = FALSE]
  }
  
  if(NROW(cluster) != n) stop("number of observations in 'cluster' and 'estfun()' do not match")

  ## catch NAs in cluster -> need to be addressed in the model object by the user
  if(anyNA(cluster)) stop("cannot handle NAs in 'cluster': either refit the model without the NA observations in 'cluster' or impute the NAs")

  ## for multi-way clustering: set up interaction patterns
  p <- NCOL(cluster)
  if (p > 1L) {
    cl <- lapply(1L:p, function(i) combn(1L:p, i, simplify = FALSE))
    cl <- unlist(cl, recursive = FALSE)
    sign <- sapply(cl, function(i) (-1L)^(length(i) + 1L))    
    paste_ <- function(...) paste(..., sep = "_")
    for (i in (p + 1L):length(cl)) {
      cluster <- cbind(cluster, Reduce(paste_, unclass(cluster[, cl[[i]] ]))) ## faster than: interaction()
    }
    if(multi0) cluster[[length(cl)]] <- 1L:n
  } else {
    cl <- list(1)
    sign <- 1
  }
  
  ## number of clusters (and cluster interactions)
  g <- sapply(1L:length(cl), function(i) {
    if(is.factor(cluster[[i]])) {
      length(levels(cluster[[i]]))
    } else {
      length(unique(cluster[[i]]))
    }
  })
  #gmin <- min(g[1L:p])
    ## FIXME: additional argument for optionally using only smallest number of clusters?
    ## See also Cameron, Gelbach and Miller (2011, page 241)
  #if(FALSE) g[] <- gmin

  ## type of model
  is_lm <- function(object) class(object)[1L] == "lm"
  is_glm <- function(object) class(object)[1L] == "glm"

  ## type of bias correction
  if(is.null(type)) {
    type <- if(is_lm(x)) "HC1" else "HC0"
  }
  type <- match.arg(type, c("HC", "HC0", "HC1", "HC2", "HC3"))
  if(type == "HC") type <- "HC0"

  ## building blocks for HC2/HC3
  if(type %in% c("HC2", "HC3"))
  {
    if(any(g == n)) h <- hatvalues(x)

    if(!all(g == n)) {
      if(!(is_lm(x) || is_glm(x))) warning("clustered HC2/HC3 are only applicable to (generalized) linear regression models")

      ## regressor matrix
      X <- model.matrix(x)
      if(any(alias <- is.na(coef(x)))) X <- X[, !alias, drop = FALSE]
      attr(X, "assign") <- NULL

      ## working weights
      w <- weights(x, "working")

      ## (X'X)^(-1)
      XX1 <- if(is.null(w)) chol2inv(qr.R(qr(X))) else chol2inv(qr.R(qr(X * sqrt(w))))

      ## working residuals
      res <- rowMeans(ef/X, na.rm = TRUE)
      res[apply(abs(ef) < .Machine$double.eps, 1L, all)] <- 0
    }
  }
  
  ## add OPG for each cluster-aggregated estfun
  for (i in 1L:length(cl))
  {
    ## estimating functions for aggregation by i-th clustering variable
    efi <- ef

    ## add cluster adjustment g/(g - 1) or not?
    ## only exception: HC0 adjustment for multiway clustering at "last interaction"
    adj <- if(multi0 & (i == length(cl))) {
               if(type == "HC1") (n - k)/(n - 1L) else 1
           } else {
               if(cadjust) g[i]/(g[i] - 1L) else 1
           }
      
    ## HC2/HC3
    if(type %in% c("HC2", "HC3")) {
      if(g[i] == n) {
        efi <- if(type == "HC2") {
	  efi/sqrt(1 - h)
	} else {
	  efi/(1 - hatvalues(x))
	}
      } else {
        for(j in unique(cluster[[i]])) {
	  ij <- which(cluster[[i]] == j)
	  Hij <- if(is.null(w)) {
	    X[ij, , drop = FALSE] %*% XX1 %*% t(X[ij, , drop = FALSE])
	  } else {
	    X[ij, , drop = FALSE] %*% XX1 %*% t(X[ij, , drop = FALSE]) %*% diag(w[ij], nrow = length(ij), ncol = length(ij))
	  }
          Hij <- if(type == "HC2") {
	    matrixpower(diag(length(ij)) - Hij, -0.5)
	  } else {
	    solve(diag(length(ij)) - Hij)
	  }
          efi[ij, ] <- drop(Hij %*% res[ij]) * X[ij, , drop = FALSE]
        }
      }

      ## "inverse" cluster adjustment that Bell & McCaffrey (2002) and hence also
      ## Cameron & Miller (2005, Eq. 25) recommend for HC3 (but not HC2)
      ## -> canceled out again if cadjust = TRUE
      efi <- sqrt((g[i] - 1L)/g[i]) * efi
    }

    ## aggregate within cluster levels      
    efi <- if(g[i] < n) apply(efi, 2L, rowsum, cluster[[i]]) else efi

    ## aggregate across cluster variables
    rval <- rval + sign[i] * adj * crossprod(efi)/n
  }

  ## HC1 adjustment with residual degrees of freedom: (n - 1)/(n - k)
  if(type == "HC1") rval <- (n - 1L)/(n - k) * rval

  return(rval)
}

## matrix power (for square root and inverse square root)
matrixpower <- function(X, p, symmetric = NULL, tol = .Machine$double.eps^(1/1.3)) {
  if((ncol(X) == 1L) && (nrow(X) == 1L)) return(X^p)
  if(is.null(symmetric)) symmetric <- isSymmetric(X)
  Xeig <- eigen(X, symmetric = symmetric)
  if(is.complex(Xeig$values)) {
    if(any(abs(Im(Xeig$values)) > tol)) warning("complex eigen values of X")
    Xeig$values <- Re(Xeig$values)
    Xeig$vectors <- Re(Xeig$vectors)
  }
  Xeig$values[Xeig$values < tol] <- 0
  # if(any(Xeig$values < 0)) stop("matrix is not positive semidefinite")
  if(symmetric) {
    Xeig$vectors %*% ((Xeig$values^p) * t(Xeig$vectors))
  } else {
    Xeig$vectors %*% ((Xeig$values^p) * matrixinverse(Xeig$vectors))
  }
}

matrixinverse <- function(X, tol = .Machine$double.eps^(1/1.3)) {
  if((ncol(X) == 1L) && (nrow(X) == 1L)) return(1/X)
  inv <- try(solve(X), silent = TRUE)
  if(!inherits(inv, "try-error")) return(inv)
  Xsvd <- svd(X)
  ok <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  inv <- Xsvd$v[, ok, drop = FALSE] %*% ((1/Xsvd$d[ok]) * t(Xsvd$u[, ok, drop = FALSE]))
  return(inv)
}
