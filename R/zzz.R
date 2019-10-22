#' nanopoRe
#'
#' methods for supporting the Oxford Nanopore Technologies Bioinformatics
#' Tutorials
#'
#' @docType package
#' @author Stephen Rudd
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @useDynLib nanopoRe, .registration = TRUE
#' @name nanopoRe
NULL

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(paste0("nanopoRe package - v ",
        utils::packageVersion("nanopoRe")))
    invisible()
}
