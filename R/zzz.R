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


.onLoad <- function(libname, pkgname) {

}

.onAttach <- function(libname, pkgname) {
    Sys.setenv("_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_"=0)
    packageStartupMessage(paste0("nanopoRe package - v ",
        utils::packageVersion("nanopoRe")))
    invisible()
}
