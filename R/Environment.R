
#' Initialise the NanopoRe environment
#'
#' Creates a NanopoRe environment; package specific parameters and values will be stored within this
#' environment; the name of the environment is defined internally
#'
#' @return None
#'
#' @examples
#' Nanopore::init()
#'
#' @export
init <- function() {
  eval(parse(text=paste0(nanopoRe.env.name," <<- new.env(parent=emptyenv())")))
  setRpath(file.path("Analysis", "R"))
}

#' check NanopoRe environment
#'
#' performs a sanity check to ensure that NanopoRe environment is initialised
#'
#' @return a logical defining whether the environment is initialised
isInitialised <- function() {
  eval(parse(text=paste0("return(exists(\"",nanopoRe.env.name, "\", mode=\"environment\"))")))
}
# nanopoRe.env.name <- "nanopoRe.env"
# usethis::use_data(nanopoRe.env.name, internal=TRUE)


getEnvironment <- function() {
  return(nanopoRe.env.name)
}
