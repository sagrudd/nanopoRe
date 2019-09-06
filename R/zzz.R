#' @useDynLib nanopoRe, .registration = TRUE
#' @importFrom Rcpp sourceCpp


.onLoad <- function(libname, pkgname) {

}

.onAttach <- function(libname, pkgname) {
    Sys.setenv("_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_"=0)
    packageStartupMessage(paste0("nanopoRe package - v ", utils::packageVersion("nanopoRe")))
    invisible()
}
