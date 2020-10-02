vcovPC <- function(x, cluster = NULL, order.by = NULL, pairwise = FALSE, sandwich = TRUE, fix = FALSE, ...)
{
  ## compute meat of sandwich
  rval <- meatPC(x, cluster = cluster, order.by = order.by, pairwise = pairwise, ...)
    
  ## full sandwich
  if (sandwich) rval <- sandwich(x, meat. = rval)

  ## check (and fix) if sandwich is not positive semi-definite
  if (fix && any((eig <- eigen(rval, symmetric = TRUE))$values < 0)) {
    eig$values <- pmax(eig$values, 0)
    rval[] <- crossprod(sqrt(eig$values) * t(eig$vectors))
  }
  return(rval)
}

meatPC <- function(x, cluster = NULL, order.by = NULL, pairwise = FALSE, kronecker = TRUE, ...)
{
  ## extract estimating functions / aka scores
  if (is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"
  ef <- estfun(x, ...)
  n <- NROW(ef)

  ## FIXME: The return value only has to be set up this way if 
  ## it is filled componentwise. But at the moment the whole matrix is
  ## computed in one big matrix computation.
  ## rval <- matrix(0, nrow = NCOL(ef), ncol = NCOL(ef),
  ##   dimnames = list(colnames(ef), colnames(ef)))

  ## cluster can either be supplied explicitly or
  ## be an attribute of the model  
  if(is.null(cluster)) cluster <- attr(x, "cluster")

  ## cluster (and potentially order.by) may be a formula
  if(inherits(cluster, "formula")) {
    ## merge if both cluster and order.by are formula: ~ cluster + order.by
    if(inherits(order.by, "formula")) {
      cluster_orderby <- ~ cluster + order.by
      cluster_orderby[[2L]][[3L]] <- order.by[[2L]]
      cluster_orderby[[2L]][[2L]] <- cluster[[2L]]
      cluster <- cluster_orderby
      order.by <- NULL
    }
    ## get variable(s) from expanded model frame
    cluster_tmp <- if("Formula" %in% loadedNamespaces()) { ## FIXME to suppress potential warnings due to | in Formula
      suppressWarnings(expand.model.frame(x, cluster, na.expand = FALSE))
    } else {
      expand.model.frame(x, cluster, na.expand = FALSE)
    }
    cluster <- model.frame(cluster, cluster_tmp, na.action = na.pass)

    ## handle omitted or excluded observations
    if((n != NROW(cluster)) && !is.null(x$na.action) && (class(x$na.action) %in% c("exclude", "omit"))) {
      cluster <- cluster[-x$na.action, , drop = FALSE]
    }
  }

  ## cluster can also be a list with both indexes
  if(is.list(cluster)) {
    if(length(cluster) > 1L & is.null(order.by)) order.by <- cluster[[2L]]
    cluster <- cluster[[1L]]
  }

  ## longitudinal time variable
  if(is.null(order.by)) order.by <- attr(x, "order.by")
  if(is.null(order.by)) {
    ix <- make.unique(as.character(cluster))
    ix <- strsplit(ix, ".", fixed = TRUE)
    ix <- sapply(ix, "[", 2L)
    ix[is.na(ix)] <- "0"
    order.by <- as.integer(ix) + 1L
  }
     
  ## handle omitted or excluded observations in both cluster and/or order.by
  ## (again, in case the formula interface was not used)
  if(any(n != c(NROW(cluster), NROW(order.by))) && !is.null(x$na.action) && class(x$na.action) %in% c("exclude", "omit")) {
    if(n != NROW(cluster))   cluster <-  cluster[-x$na.action]
    if(n != NROW(order.by)) order.by <- order.by[-x$na.action]
  }

  ## model matrix
  X <- model.matrix(x)
 
  if (any(alias <- is.na(coef(x)))) X <- X[, !alias, drop = FALSE]
  attr(X, "assign") <- NULL
    
  ## working residuals
  res <- rowMeans(ef/X, na.rm = TRUE)
  res[apply(abs(ef) < .Machine$double.eps, 1L, all)] <- 0

    ## balanced panel
    balanced <- isTRUE(length(unique(table(cluster, order.by))) == 1L)
    if(balanced) {

    ## split residuals by cluster
    e <- split(res, cluster)    

    ## bind into matrix
    e <- do.call("cbind", e)
        
    ## set up omega
    N <- prod(dim(e))
    t <- dim(e)[1L]
    n <- dim(e)[2L]
    sigma <- crossprod(e) / t
    if(kronecker) {
        omega <- kronecker(sigma, diag(1L, t))
        omega <- t(X) %*% omega %*% X / N
    } else {
        xx <- split(X, cluster)
        xxCL <- lapply(1L:n, function(i) matrix(xx[[i]], t, ncol(X)))
        omega <- Reduce("+", lapply(1L:n, function(j) Reduce("+", lapply(1L:n, function(i) sigma[j,i] * t(xxCL[[j]]) %*% xxCL[[i]] / N))))
    }
    }

    
    ## unbalanced panel  
    if(!balanced) {
        pair <- data.frame(cluster = cluster, order.by = order.by)
        num <- dim(pair)[1L]
        pair$res <- res
        pair <- data.frame(pair, X)
        Tij <- expand.grid(cluster = unique(cluster), order.by = unique(order.by))
        pair <- merge(pair, Tij, by = c("cluster", "order.by"), all = TRUE)
        
    ## balance panel
       if(!pairwise) {
           
    ## extract "full" clusters   
    rem <- subset(pair$order.by, is.na(pair$res))
    pairX <- as.matrix(pair[, 4L:dim(pair)[2L]])    
    pairX[is.na(pairX)] <- 0
    pair[which(pair$order.by %in% rem), 3L:dim(X)[2L]] <- NA
    
    pairNA <- na.omit(pair)
    pair[is.na(pair)] <- 0
    res <- pairNA[, 3L]
    clpairNA <- pairNA[, 1L]
    clpair <- pair[, 1L]
           
    ## split residuals by cluster
    e <- split(res, clpairNA)    

    ## bind into matrix
    e <- do.call("cbind", e)
        
    ## set up omega
    t <- dim(e)[1L]
    n <- dim(e)[2L]
    tt <- length(unique(order.by))
    sigma <- crossprod(e) / t

    if(kronecker) {
        omega <- kronecker(sigma, diag(1L, tt))
        omega <- t(pairX) %*% omega %*% pairX / num
    } else {
        xx <- split(pairX, clpair)
        xxCL <- lapply(1L:n, function(i) matrix(xx[[i]], tt, ncol(X)))
        omega <- Reduce("+", lapply(1L:n, function(j) Reduce("+", lapply(1L:n, function(i) sigma[j,i] * t(xxCL[[j]]) %*% xxCL[[i]] / num))))
    }
    } else {
        
    ## use pairwise calculation for omega (Bailey and Katz, 2011)
    N <- dim(X)[1L]    
    e <- pair[, 3L]
    X <- pair[, 4L:dim(pair)[2L]]
   
    narows <- apply(!is.na.data.frame(X), 1, prod)
    namatrix <- matrix(narows, length(unique(cluster)), length(unique(order.by)), byrow = TRUE)
    namatrix <- t(namatrix)
    denoma <- crossprod(namatrix)
    X[is.na(X)] <- 0L
    
    ## set up omega
    e <- matrix(e, length(unique(cluster)), length(unique(order.by)), byrow = TRUE)
    e[is.na(e)] <- 0L
    e <- t(e)
    sigma <- crossprod(e) / denoma
    if(kronecker) {
        omega <- kronecker(sigma, diag(1L, NROW(e)))
        X <- as.matrix(X)
        X[is.na(X)] <- 0L
        omega <- t(X) %*% omega %*% X / N
    } else {
        t <- length(unique(pair$order.by))
        n <- length(unique(pair$cluster))
        xx <- split(X, pair$cluster)
        xxCL <- lapply(1L:n, function(i) as.matrix(xx[[i]], t, ncol(X)))
        omega <- Reduce("+", lapply(1L:n, function(j) Reduce("+", lapply(1L:n, function(i) sigma[j,i] * t(xxCL[[j]]) %*% xxCL[[i]] / N))))
    }
    }
    }
    
    rval <- omega
    
  return(rval = rval)
}
