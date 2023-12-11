vcovBS.lm <- function(x, cluster = NULL, R = 250, type = "xy", ..., fix = FALSE, use = "pairwise.complete.obs", applyfun = NULL, cores = NULL, qrjoint = FALSE, center = "mean")
{
  ## set up return value with correct dimension and names
  cf0 <- coef(x)
  k <- length(cf0)
  n <- nobs(x)
  rval <- matrix(0, nrow = k, ncol = k, dimnames = list(names(cf0), names(cf0)))

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
  
  if(NROW(cluster) != n) stop("number of observations in 'cluster' and 'nobs()' do not match")

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
  } else {
    cl <- list(1)
    sign <- 1
  }

  ## bootstrap type: xy vs. residual vs. wild (in various flavors)
  if(is.function(type)) {
    wild <- type
    type <- "user"
  } else {
    wild <- function(n) NULL
  }
  type <- match.arg(gsub("wild-", "", tolower(type), fixed = TRUE),
    c("xy", "jackknife", "fractional", "residual", "wild", "webb", "rademacher", "mammen", "norm", "user"))
  if(type == "wild") type <- "rademacher"
  
  ## set up wild bootstrap function
  wild <- switch(type,
    "rademacher" = function(n) sample(c(-1, 1), n, replace = TRUE),
    "mammen"     = function(n) sample(c(-1, 1) * (sqrt(5) + c(-1, 1))/2, n, replace = TRUE,
				 prob = (sqrt(5) + c(1, -1))/(2 * sqrt(5))),
    "norm"       = function(n) rnorm(n),
    "webb"       = function(n) sample(c(-sqrt((3:1)/2), sqrt((1:3)/2)), n, replace = TRUE),
    wild
  )

  ## model information: original response, weights, and design matrix or the corresponding fitted/residuals/QR
  y <- if(!is.null(x$y)) {
    x$y
  } else if(!is.null(x$model)) {
    model.response(x$model)
  } else {
    model.response(model.frame(x))
  }
  wts <- if(is.null(x$weights)) rep.int(1, n) else x$weights
  off <- if(is.null(x$offset))  rep.int(0, n) else x$offset
  xfit <- if(type %in% c("xy", "jackknife", "fractional")) model.matrix(x) else list(fit = x$fitted.values, res = x$residuals, qr = x$qr)

  ## apply infrastructure for refitting models
  if(is.null(applyfun)) {
    applyfun <- if(is.null(cores)) {
      lapply
    } else {
      if(.Platform$OS.type == "windows") {
        cl_cores <- parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl_cores))
        function(X, FUN, ...) parallel::parLapply(cl = cl_cores, X, FUN, ...)
      } else {
        function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = cores)
      }
    }
  }

  ## bootstrap for each cluster dimension
  for (i in 1L:length(cl)) {
    ## cluster structure
    cli <- if(type %in% c("xy", "jackknife", "residual")) {
      split(seq_along(cluster[[i]]), cluster[[i]])
    } else {
      factor(cluster[[i]], levels = unique(cluster[[i]]))
    }

    ## check for balanced clusters in residual bootstrap
    if(type == "residual" && length(unique(sapply(cli, length))) > 1L) {
      warning("residual bootstrap is not well-defined for unbalanced clusters")
    } 

    ## bootstrap fitting function
    bootfit <- switch(type,
      "xy" = function(j, ...) {
        j <- unlist(cli[sample(names(cli), length(cli), replace = TRUE)])
        .lm.fit(xfit[j, , drop = FALSE] * sqrt(wts[j]), (y[j] - off[j]) * sqrt(wts[j]), ...)$coefficients
      },
      "jackknife" = function(j, ...) {
        j <- unlist(cli[-j])
        .lm.fit(xfit[j, , drop = FALSE] * sqrt(wts[j]), (y[j] - off[j]) * sqrt(wts[j]), ...)$coefficients
      },
      "fractional" = function(j, ...) {
        fw <- rexp(nlevels(cli))
        fw <- fw[cli]/mean(fw)
        .lm.fit(xfit * sqrt(wts * fw), (y - off) * sqrt(wts * fw), ...)$coefficients
      },
      "residual" = function(j, ...) {
        j <- unlist(cli[sample(names(cli), length(cli), replace = TRUE)])
        yboot <- xfit$fit + xfit$res[j]
        yboot <- (yboot - off) * sqrt(wts)
	if(qrjoint) yboot else qr.coef(xfit$qr, yboot)
      },
      function(j, ...) {
        j <- wild(nlevels(cli))
        yboot <- xfit$fit + xfit$res * j[cli]
        yboot <- (yboot - off) * sqrt(wts)
	if(qrjoint) yboot else qr.coef(xfit$qr, yboot)
      }
    )

    ## for jackknife the number of replications is always the number of "observations" (cluster units)
    if(type == "jackknife") R <- length(cli)
    
    ## actually refit
    cf <- applyfun(1L:R, bootfit, ...)

    ## aggregate across cluster variables
    if(type == "jackknife") {
      cf <- do.call("cbind", cf)
      center <- match.arg(center, c("mean", "estimate"))
      rval <- rval + sign[i] * (R - 1L)/R * tcrossprod(cf - if(center == "mean") rowMeans(cf) else cf0)
    } else {
      cf <- if(qrjoint && type != "xy") t(qr.coef(xfit$qr, do.call("cbind", cf))) else do.call("rbind", cf)
      rval <- rval + sign[i] * cov(cf, use = use)
    }
  }

  ## check (and fix) if sandwich is not positive semi-definite
  if(fix && any((eig <- eigen(rval, symmetric = TRUE))$values < 0)) {
    eig$values <- pmax(eig$values, 0)
    rval[] <- crossprod(sqrt(eig$values) * t(eig$vectors))
  }
  return(rval)
}

vcovBS.glm <- function(x, cluster = NULL, R = 250, start = FALSE, type = "xy", ..., fix = FALSE, use = "pairwise.complete.obs", applyfun = NULL, cores = NULL, center = "mean")
{
  ## set up return value with correct dimension and names
  cf0 <- coef(x)
  if(identical(start, TRUE))  start <- cf0
  if(identical(start, FALSE)) start <- NULL
  k <- length(cf0)
  n <- nobs(x)
  rval <- matrix(0, nrow = k, ncol = k, dimnames = list(names(cf0), names(cf0)))

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
  
  if(NROW(cluster) != n) stop("number of observations in 'cluster' and 'nobs()' do not match")

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
  } else {
    cl <- list(1)
    sign <- 1
  }

  ## model information: original response, design matrix, fitting method
  y <- if(!is.null(x$y)) {
    x$y
  } else if(!is.null(x$model)) {
    model.response(x$model)
  } else {
    model.response(model.frame(x))
  }
  xfit <- model.matrix(x)
  method <- x$method

  ## apply infrastructure for refitting models
  if(is.null(applyfun)) {
    applyfun <- if(is.null(cores)) {
      lapply
    } else {
      if(.Platform$OS.type == "windows") {
        cl_cores <- parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl_cores))
        function(X, FUN, ...) parallel::parLapply(cl = cl_cores, X, FUN, ...)
      } else {
        function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = cores)
      }
    }
  }

  ## xy bootstrap vs. jackknife
  type <- match.arg(tolower(type), c("xy", "jackknife", "fractional"))

  ## bootstrap for each cluster dimension
  for (i in 1L:length(cl))
  {
    ## cluster structure
    cli <- if(type != "fractional") {
      split(seq_along(cluster[[i]]), cluster[[i]])
    } else {
      factor(cluster[[i]], levels = unique(cluster[[i]]))
    }

    ## bootstrap fitting function
    bootfit <- function(j, ...) {
      if(type == "xy") {
        j <- unlist(cli[sample(names(cli), length(cli), replace = TRUE)])
        wts <- 1
      } else if(type == "jackknife") {
        j <- unlist(cli[-j])
        wts <- 1
      } else if(type == "fractional") {
        j <- 1L:n
        wts <- rexp(nlevels(cli))
        wts <- wts[cli]/mean(wts)
      }
      eval(
        call(if(is.function(method)) "method" else method, 
          x = xfit[j, , drop = FALSE], y = y[j], weights = x$prior.weights[j] * wts, offset = x$offset[j],
          family = x$family, start = start, control = x$control, intercept = attr(x$terms, "intercept") > 0)
      )$coefficients
    }

    ## for jackknife the number of replications is always the number of "observations" (cluster units)
    if(type == "jackknife") R <- length(cli)
    
    ## actually refit
    cf <- applyfun(1L:R, bootfit, ...)

    ## aggregate across cluster variables
    if(type == "jackknife") {
      cf <- do.call("cbind", cf)
      center <- match.arg(center, c("mean", "estimate"))
      rval <- rval + sign[i] * (R - 1L)/R * tcrossprod(cf - if(center == "mean") rowMeans(cf) else cf0)
    } else {
      cf <- do.call("rbind", cf)
      rval <- rval + sign[i] * cov(cf, use = use)
    }
  }

  ## check (and fix) if sandwich is not positive semi-definite
  if(fix && any((eig <- eigen(rval, symmetric = TRUE))$values < 0)) {
    eig$values <- pmax(eig$values, 0)
    rval[] <- crossprod(sqrt(eig$values) * t(eig$vectors))
  }
  return(rval)
}
