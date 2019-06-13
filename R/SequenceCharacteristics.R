
#' calculate Nstatistics (e.g. N50) for given sequence collection
#'
#' Given a set of contigs, the N50 is defined as the sequence length of the shortest contig at 50% of the total genome length (https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics)
#'
#' @param len.vector is a vector of sequence lengths
#' @param n is a numeric of the value to use
#' @return numeric of the corresponding N value
#'
#' @examples
#' # to calculate e.g. N50
#' ncalc(seq(1,10), 0.5)
#'
#' @export
ncalc <- function(len.vector, n) {
  len.sorted <- rev(sort(len.vector))
  len.sorted[cumsum(len.sorted) >= sum(len.sorted)*n][1]
}


#' calculate Lstatistics (e.g. L50) for given sequence collection
#'
#' Given a set of contigs, each with its own length, the L50 count is defined as the smallest number of contigs whose length sum makes up half of genome size (https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics)
#'
#' @param filename is a qualified path to a file of interest
#' @return character representation of md5sum
#'
#' @examples
#' lcalc(seq(1, 10), 0.5)
#'
#' @export
lcalc <- function(len.vector, n) {
  len.sorted <- rev(sort(len.vector))
  which(cumsum(len.sorted) >= sum(len.sorted)*n)[1]
}
