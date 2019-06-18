#' @useDynLib nanopoRe, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL



.onLoad <- function(libname, pkgname) {
  message(paste0("nanopoRe package - v ", packageVersion("nanopoRe")))
  message(paste0("    dependency(Biostrings)"))
  suppressMessages(library(Biostrings))
  message(paste0("    dependency(dplyr)"))
  suppressMessages(library(dplyr))
  message(paste0("    dependency(Rsamtools)"))
  suppressMessages(library(Rsamtools))
  message(paste0("    dependency(GenomicAlignments)"))
  suppressMessages(library(GenomicAlignments))
  message(paste0("    dependency(ShortRead)"))
  suppressMessages(library(ShortRead))
  init()
  invisible()
}
