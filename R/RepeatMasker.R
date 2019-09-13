
# http://www.repeatmasker.org/species/hg.html hg19 - Feb 2009 - RepeatMasker
# open-4.0.5 - Repeat Library 20140131
# http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/...
#    hg19.fa.out.gz



#' extract GRanges object corresponding to RepeatMasker annotation
#'
#' this method is used to extract most salient data from the RepeatMasker
#' annotation file into a GRanges
#' object so that Genome Geometry analyses can be performed and genome
#' coordinates can be used to quickly
#' associate with repeat types
#'
#' @usage RepeatMaskerGR(RepeatMaskerFile = NULL, RepeatMaskerSRC = NULL,
#'     force = FALSE)
#' @importFrom RCurl getBinaryURL
#' @importFrom data.table fread
#' @param RepeatMaskerFile is a path to an already downloaded RepeatMasker
#' results file
#' @param RepeatMaskerSRC is the URL to RepeatMasker dataset
#' @param force is a boolean to specify whether the build of external data
#' should be forced
#' @return GRanges object of annotated repeats
#'
#' @examples
#' # Using C.elegans repeats as demo; small file with only 90k repeats = easy
#' src <- paste0("http://www.repeatmasker.org/genomes/ce10/",
#'     "RepeatMasker-rm405-db20140131/ce10.fa.out.gz")
#' # if you roll your eyes and sigh - then please feel my frustration with the
#' # required strict adherence to Bioc conventions for line length ...
#' repeats <- RepeatMaskerGR(RepeatMaskerSRC=src)
#'
#' @export
RepeatMaskerGR <- function(RepeatMaskerFile = NULL, RepeatMaskerSRC = NULL,
    force = FALSE) {
    if (is.null(RepeatMaskerSRC)) {
        RepeatMaskerSRC = paste0("http://www.repeatmasker.org/genomes/hg19/",
        "RepeatMasker-rm405-db20140131/hg19.fa.out.gz")
    }
    if (is.null(RepeatMaskerFile)) {
        RepeatMaskerFile <- file.path(getRpath(), basename(RepeatMaskerSRC))
        if (!file.exists(RepeatMaskerFile) || force) {
            content = getBinaryURL(RepeatMaskerSRC)
            writeBin(content, RepeatMaskerFile)
        }
    }

    RepeatMaskerGRData <- file.path(getRpath(), paste0(sub("\\.[^.]*$", "",
        basename(RepeatMaskerFile)), ".GRanges", ".Rdata"))
    if (file.exists(RepeatMaskerGRData) & !force) {
        return(readRDS(file = RepeatMaskerGRData))
    }

    repeatMaskData <- data.table::fread(RepeatMaskerFile, stringsAsFactors =
        FALSE, fill = TRUE, skip = 2)
    GRdata <- GRanges(seqnames = gsub("chr", "", unlist(repeatMaskData[, 5])),
        ranges = IRanges(start = unlist(repeatMaskData[,6]), end =
        unlist(repeatMaskData[, 7])), strand = gsub("C", "-",
        unlist(repeatMaskData[, 9])),
        repeatClass = unlist(repeatMaskData[, 11]),
        repeatUnit = unlist(repeatMaskData[, 10]))
    saveRDS(GRdata, file = RepeatMaskerGRData)
    return(GRdata)

}
