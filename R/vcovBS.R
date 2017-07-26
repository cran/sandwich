.vcovBSenv <- new.env()

vcovBS <- function(x, ...) {
  UseMethod("vcovBS")
}

vcovBS.default <- function(x, cluster = NULL, R = 250, start = FALSE, ...)
{
  ## cluster variable (default 1:n)
  if(is.null(cluster)) {
    cluster <- try(1L:NROW(estfun(x)), silent = TRUE)
    if(inherits(cluster, "try-error")) cluster <- 1L:NROW(model.frame(x))
  }
  cl <- split(seq_along(cluster), cluster)
  
  ## set up coefficient matrix
  cf <- coef(x)
  cf <- matrix(rep.int(0, length(cf) * R), ncol = length(cf),
    dimnames = list(NULL, names(cf)))

  ## use starting values?
  assign(".vcovBSstart", if(isTRUE(start)) coef(x) else NULL, envir = .vcovBSenv)

  ## update on bootstrap samples
  for(i in 1:R) {
      subset <- unlist(cl[sample(names(cl), length(cl), replace = TRUE)])
      assign(".vcovBSsubset", subset, envir = .vcovBSenv)
      up <- if(is.null(.vcovBSenv$.vcovBSstart)) {
        update(x, subset = .vcovBSenv$.vcovBSsubset, ..., evaluate = FALSE)
      } else {
        update(x, subset = .vcovBSenv$.vcovBSsubset, start = .vcovBSenv$.vcovBSstart, ..., evaluate = FALSE)      
      }
      env <- try(environment(terms(x)))
      if(inherits(env, "try-error")) env <- NULL
      up <- eval(up, envir = env, enclos = parent.frame())
      remove(".vcovBSsubset", envir = .vcovBSenv)
      cf[i, ] <- coef(up)
  }
  remove(".vcovBSstart", envir = .vcovBSenv)
  
  # covariance of bootstrap coefficients
  cov(cf)
}

