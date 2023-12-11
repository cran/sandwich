.vcovBSenv <- new.env()

vcovBS <- function(x, ...) {
  UseMethod("vcovBS")
}

vcovBS.default <- function(x, cluster = NULL, R = 250, start = FALSE, type = "xy", ..., fix = FALSE, use = "pairwise.complete.obs", applyfun = NULL, cores = NULL, center = "mean")
{
  ## set up return value with correct dimension and names
  cf0 <- coef(x)
  k <- length(cf0)
  n <- nobs0(x)
  rval <- matrix(0, nrow = k, ncol = k, dimnames = list(names(cf0), names(cf0)))

  ## try to figure out environment for update()
  env <- try(environment(terms(x)))
  if(inherits(env, "try-error")) env <- NULL

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

  ## use starting values?
  assign(".vcovBSstart", if(isTRUE(start)) coef(x) else NULL, envir = .vcovBSenv)

  ## xy bootstrap vs. jackknife
  type <- match.arg(tolower(type), c("xy", "jackknife", "fractional"))

  ## bootstrap for each cluster dimension
  for (i in 1L:length(cl)) {
    ## cluster structure
    cli <- if(type != "fractional") {
      split(seq_along(cluster[[i]]), cluster[[i]])    
    } else {
      factor(cluster[[i]], levels = unique(cluster[[i]]))
    }

    ## bootstrap fitting function via update()
    bootfit <- function(j, ...) {
        clj <- switch(type,
          "xy" = {
            unlist(cli[sample(names(cli), length(cli), replace = TRUE)])
          },
          "jackknife" = {
            unlist(cli[-j])
          },
          "fractional" = {
            fw <- rexp(nlevels(cli))
            fw <- fw[cli]/mean(fw)
            mw <- weights(x)
            j <- if(is.null(mw)) fw else mw * fw
          })
        assign(".vcovBSsubset", clj, envir = .vcovBSenv)
        up <- if(is.null(.vcovBSenv$.vcovBSstart)) {
          if(type != "fractional") {
            update(x, subset = .vcovBSenv$.vcovBSsubset, ..., evaluate = FALSE)
          } else {
            update(x, weights = .vcovBSenv$.vcovBSsubset, ..., evaluate = FALSE)
          }
        } else {
          if(type != "fractional") {
            update(x, subset = .vcovBSenv$.vcovBSsubset, start = .vcovBSenv$.vcovBSstart, ..., evaluate = FALSE)      
          } else {
            update(x, weights = .vcovBSenv$.vcovBSsubset, start = .vcovBSenv$.vcovBSstart, ..., evaluate = FALSE)                
          }
        }
        up <- eval(up, envir = env, enclos = parent.frame())
        coef(up)
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
  ## clean up starting values again
  remove(".vcovBSstart", envir = .vcovBSenv)

  ## check (and fix) if sandwich is not positive semi-definite
  if(fix && any((eig <- eigen(rval, symmetric = TRUE))$values < 0)) {
    eig$values <- pmax(eig$values, 0)
    rval[] <- crossprod(sqrt(eig$values) * t(eig$vectors))
  }
  return(rval)
}
