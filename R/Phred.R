
#' calculate mean Phred scores from list of Q values
#'
#' This method calculates the mean quality value from a list of q-values;
#' handles the challenge of relative log-scaling of the Q and challenges of a
#' more linear mean; this avoids an artificial inflation of phred-scores
#'
#' @param q is a vector of q values
#' @return mean phred scaled q-value
#'
#' @examples
#' mean(c(20,30,40))
#' phredmean(c(20,30,40))
#'
#' @export
phredmean <- function(q) {
    -10 * log10(mean(10^(q/-10), na.rm = TRUE))
}



#' calculate mean Phred score from an ASCII encoded phred string
#'
#' FASTQ and BAM store per base qualities as an ASCII string. Accessory methods
#' in e.g. ShortRead allow for a sum of the numeric encoded scores; this is not
#' corrected for the log/linear so scores are synthetically boosted
#' - this simple method performs mean on the character level data ...
#'
#' @param qstr is an ASCII encoded Phred quality score
#' @return mean phred scaled q-value
#'
#' @examples
#' qualToMeanQ('ABCDEF')
#'
#' @export
qualToMeanQ <- function(qstr) {
    baseq <- as.numeric(charToRaw(qstr)) - 33
    meanerror <- mean(10^(baseq/-10))
    -10 * log10(meanerror)
}
