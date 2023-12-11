vcovJK <- function(x, ...) {
  UseMethod("vcovJK")
}

vcovJK.default <- function(x, cluster = NULL, center = "mean", ...) {
  center <- match.arg(center, c("mean", "estimate"))
  vcovBS(x, cluster = cluster, center = center, type = "jackknife", ...)
}
