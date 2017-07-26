vcovPL <- function(x, cluster = NULL, order.by = NULL,
  kernel = "Bartlett", sandwich = TRUE, fix = FALSE, ...)
{
  ## compute meat of sandwich
  rval <- meatPL(x, cluster = cluster, order.by = order.by, kernel = kernel, ...)

  ## full sandwich
  if(sandwich) rval <- sandwich(x, meat. = rval)

  ## check (and fix) if sandwich is not positive semi-definite
  if(fix && any((eig <- eigen(rval, symmetric = TRUE))$values < 0)) {
    eig$values <- pmax(eig$values, 0)
    rval <- crossprod(sqrt(eig$values) * t(eig$vectors))
  }
  return(rval)
}

meatPL <- function(x, cluster = NULL, order.by = NULL,
  kernel = "Bartlett", lag = "NW1987", bw = NULL, adjust = TRUE, ...) ## adjust/cadjust?
{
  ## extract estimating functions / aka scores
  if (is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"
  ef <- estfun(x, ...)
  k <- NCOL(ef)
  n <- NROW(ef)

  ## cluster can either be supplied explicitly or
  ## be an attribute of the model...FIXME: other specifications?
  if(is.null(cluster)) cluster <- attr(x, "cluster")

  ## cluster can also be a list with both indexes
  if(is.list(cluster)) {
    if(length(cluster) > 1L & is.null(order.by)) order.by <- cluster[[2L]]
    cluster <- cluster[[1L]]
  }

  ## Newey-West type standard errors if neither cluster nor order.by is specified
  if(is.null(cluster) && is.null(order.by)) cluster <- rep(1, n)
        
  ## resort to cross-section if no clusters are supplied
  if(is.null(cluster)) cluster <- 1L:n
    
  ## longitudinal time variable
  if(is.null(order.by)) order.by <- attr(x, "order.by")
  if(is.null(order.by)) {
    ix <- make.unique(as.character(cluster))
    ix <- strsplit(ix, ".", fixed = TRUE)
    ix <- sapply(ix, "[", 2L)
    ix[is.na(ix)] <- "0"
    order.by <- as.integer(ix) + 1L
  }

  ## reorder scores and cluster/time variables
  index <- order(order.by)
  ef <- ef[index, , drop = FALSE]
  order.by <- order.by[index]
  cluster <- cluster[index]

  ## aggregate within time periods
  if(length(unique(order.by)) < n) ef <- apply(ef, 2L, tapply, order.by, sum)
  nt <- NROW(ef)
  
  ## lag/bandwidth selection
    if(is.character(lag)) {
        if(lag == "P2009") lag <- "max" 
        switch(match.arg(lag, c("max", "NW1987", "NW1994")),
               "max" = {
                   lag <- nt - 1L
               },
               "NW1987" = {
                   lag <- floor(nt^(1/4))
                       },
               "NW1994" = {
                   lag <- floor(4L * (nt / 100L)^(2/9))
                       }
               )}
  if(is.numeric(lag)) lag <- lag
  if(is.null(lag) & is.null(bw)) lag <- nt - 1L
  if(!is.null(lag) & is.null(bw)) bw <- lag + 1

  ## set up kernel weights up to maximal number of lags
  weights <- kweights(0L:(nt - 1L)/bw, kernel = kernel)

  rval <- 0.5 * crossprod(ef) * weights[1L]
    
  if(length(weights) > 1L) {
    for (ii in 2L:length(weights)) {
      rval <- rval + weights[ii] * crossprod(ef[1L:(nt - ii + 1), , drop = FALSE], ef[ii:nt, , drop = FALSE])
    }
 }
    
  rval <- rval + t(rval)

  if(adjust) rval <- n/(n - k) * rval

  return(rval/n)
}



