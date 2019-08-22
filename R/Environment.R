
#' Initialise the NanopoRe environment
#'
#' Creates a NanopoRe environment; package specific parameters and values will be stored within this
#' environment; the name of the environment is defined internally
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' Nanopore::init()
#' }
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
  if (!isInitialised()) {
    init()
  }
  return(nanopoRe.env.name)
}

getCachedFileObject <- function(objectName, fileName) {
  if (hasCachedObject(objectName)) {
    return(getCachedObject(objectName))
  } else {
    setCachedObject(objectName, readRDS(file=fileName))
    return(getCachedObject(objectName))
  }
}

setCachedObject <- function(objectName, data) {
  #message(paste0("caching object ",objectName,"\n"))
  assign(objectName, data, envir=get(getEnvironment()))
}

hasCachedObject <- function(objectName) {
  if (exists(objectName, envir=get(getEnvironment()))) {
    return(TRUE)
  }
  return(FALSE)
}

getCachedObject <- function(objectName) {
  #message(paste0("uncaching object ",objectName,"\n"))
  return(get(objectName, envir=get(getEnvironment())))
}


#' Initialise the NanopoRe environment
#'
#' Creates a NanopoRe environment; package specific parameters and values will be stored within this
#' environment; the name of the environment is defined internally
#'
#' @importFrom utils ls.str
#' @return None
#'
#' @examples
#' \dontrun{
#' nanopoRe::init()
#' }
#' @export
listCachedObjects <- function() {
  ls.str(get(getEnvironment()))
}
