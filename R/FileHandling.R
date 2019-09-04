

#' perform a sanity check and prepare summary info on fastq file
#'
#' This method will prepare summary statistics on a fastq file and cache results (assuming immutable nature)
#'
#' @importFrom gdata humanReadable
#' @importFrom scales comma_format
#' @param input_fastq is the path to the fastq file to check
#' @param force logical value describing whether the analysis should be force recalculated
#' @return vector of observations
#'
#' @examples
#' fastqFile <- system.file("extdata", "example.fastq.gz", package = "nanopoRe")
#' fqcheck <- fastqCheckup(fastqFile)
#' names(fqcheck)
#' fqcheck[['reads']]
#'
#' @export
fastqCheckup <- function(input_fastq, force = FALSE) {
    fileResults <- file.path(getRpath(), paste0(sub("\\.[^.]*$", "", basename(input_fastq)), ".sanityfq",
        ".Rdata"))
    if (file.exists(fileResults) & !force) {
        return(readRDS(file = fileResults))
    }
    if (!fastqValidator(input_fastq)) {
        # something is FUBAR with this file ...
        saveRDS(FALSE, file = fileResults)
        return(FALSE)
    }

    bases <- gdata::humanReadable(getFastqBases(), standard = "SI")
    bases <- gsub("GB", "Gbases", bases)
    bases <- gsub("MB", "Mbases", bases)
    bases <- gsub("kB", "kbases", bases)

    fileinfo <- c(filename = basename(input_fastq), checksum = md5sum(input_fastq), reads = (scales::comma_format())(getFastqCount()),
        bases = bases, size = gdata::humanReadable(file.info(input_fastq)$size), nt = getFastqBases(),
        n = getFastqCount())
    saveRDS(fileinfo, file = fileResults)
    return(fileinfo)
}
