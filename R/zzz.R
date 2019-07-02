#' @useDynLib nanopoRe, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL



.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("nanopoRe package - v ", utils::packageVersion("nanopoRe")))
  init()
  invisible()
}
