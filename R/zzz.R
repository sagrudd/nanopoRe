#' @useDynLib nanopoRe, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL



.onLoad <- function(libname, pkgname) {
  packageStartupMessage(paste0("nanopoRe package - v ", utils::packageVersion("nanopoRe")))
  init()
  invisible()
}
