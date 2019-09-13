
#' calculate md5sum checksum for qualified file
#'
#' This method calculates the md5 checksum for the specified file
#'
#' @importFrom digest digest
#' @param filename is a qualified path to a file of interest
#' @return character representation of md5sum
#'
#' @examples
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2",
#'     package = "nanopoRe")
#' md5sum(seqsumFile)
#'
#' @export
md5sum <- function(filename) {
    digest(filename, algo = "md5", file = TRUE)
}
