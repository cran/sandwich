.First.lib <- function(lib, pkg) {
  if(!("package:zoo" %in% search() || require(zoo))) warning("could not load package zoo")
}
