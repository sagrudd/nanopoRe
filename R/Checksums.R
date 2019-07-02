
#' calculate md5sum checksum for qualified file
#'
#' This method calculates the md5 checksum for the specified file
#'
#' @importFrom digest digest
#' @param filename is a qualified path to a file of interest
#' @return character representation of md5sum
#'
#' @examples
#' \dontrun{
#' md5sum(file.path("DESCRIPTION"))
#' }
#'
#' @export
md5sum <- function(filename) {
  digest(filename, algo="md5", file=TRUE)
}
