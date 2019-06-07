
#' set path for R results
#'
#' This method sets a path where default R results will be stored during an analysis
#'
#' @param path is the location to the directory for storing results
#' @return None
#'
#' @examples
#' setRpath(file.path("Analysis", "R"))
#'
#' @export
setRpath <- function(path) {
  assign("r_path", path, envir=get(getEnvironment()))
}



#' get path for R results
#'
#' This method returns the path where default R results will be stored during an analysis
#'
#' @return character representation of path
#'
#' @examples
#' getRpath()
#'
#' @export
getRpath <- function(path) {
  return(get("r_path", envir=get(getEnvironment())))
}
