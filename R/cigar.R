

#' prepares a matrix of operations and coordinates from a CIGAR
#'
#' This method will prepare a matrix of operations, query coordinates
#' and target coordinates from a Cigar string extracted from a BAM file
#'
#' @param cigar is the CIGAR string to parse
#'
#' @return matrix of coordinates
#'
#' @examples
#' cigarStr <-
#'     "46S14M5D14M2D44M3I6M2I14M3I15M1D27M1D11M1I20M1I8M2I2M3I10S"
#' cigarPM <- buildPositionMap(cigarStr)
#'
#' @export
buildPositionMap <- function(cigar) {
    # this method is intended only as an interim solution until
    ###### GenomicAlignments::queryLoc2refLoc has been implemented (@sagrudd - 21/10/2019)

    # BIO::Cigar from https://github.com/MullinsLab/Bio-Cigar/blob/master/lib/Bio/Cigar.pm
    # https://github.com/samtools/htslib/blob/develop/htslib/sam.h#L81-L100

    op_consumes <- data.frame(query=c(1,1,0,0,1,0,0,1,1),
                              reference=c(1,0,1,1,0,0,0,1,1),
                            row.names=c('M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'))

    rle <- GenomicAlignments::cigarToRleList(cigar)[[1]]
    rleEvents <- rep(rle@values, rle@lengths)
    qpos <- cumsum(op_consumes[rleEvents,'query'])
    rpos <- cumsum(op_consumes[rleEvents,'reference'])
    positionMap <- data.frame(event = rleEvents,
                              qpos = qpos,
                              rpos = rpos)
    return(positionMap)
}


#' maps reference coordinates back to the query sequence
#'
#' This method will prepare a vector of query coordinates from
#' the supplied reference coordinates
#'
#' @param poi is a vector of numeric positions of interest
#' @param cigar is the CIGAR string to parse
#'
#' @return matrix of coordinates
#'
#' @examples
#' cigarStr <-
#'     "46S14M5D14M2D44M3I6M2I14M3I15M1D27M1D11M1I20M1I8M2I2M3I10S"
#' cigarR2Q(c(20,40,60,80,100), cigarStr)
#'
#' @export
cigarR2Q <- function(poi, cigar) {
    pmap <- buildPositionMap(cigar)
    return(pmap[match(poi, pmap$rpos),"qpos"])
}


#' maps query coordinates to the reference sequence
#'
#' This method will prepare a vector of reference coordinates from
#' the supplied query coordinates
#'
#' @param poi is a vector of numeric positions of interest
#' @param cigar is the CIGAR string to parse
#'
#' @return matrix of coordinates
#'
#' @examples
#' cigarStr <-
#'     "46S14M5D14M2D44M3I6M2I14M3I15M1D27M1D11M1I20M1I8M2I2M3I10S"
#' cigarQ2R(c(20,40,60,80,100), cigarStr)
#'
#' @export
cigarQ2R <- function(poi, cigar) {
    pmap <- buildPositionMap(cigar)
    return(pmap[match(poi, pmap$qpos),"rpos"])
}

#' extracts cigar operations based on reference coordinates
#'
#' This method will prepare a factor of cigar operations from
#' the supplied reference coordinates
#'
#' @param poi is a vector of numeric positions of interest
#' @param cigar is the CIGAR string to parse
#'
#' @return factor of operations
#'
#' @examples
#' cigarStr <-
#'     "46S14M5D14M2D44M3I6M2I14M3I15M1D27M1D11M1I20M1I8M2I2M3I10S"
#' cigarR2Op(c(20,40,60,80,100), cigarStr)
#'
#' @export
cigarR2Op <- function(poi, cigar) {
    pmap <- buildPositionMap(cigar)
    return(pmap[match(poi, pmap$rpos),"event"])
}

#' extracts cigar operations based on query coordinates
#'
#' This method will prepare a factor of cigar operations from
#' the supplied query coordinates
#'
#' @param poi is a vector of numeric positions of interest
#' @param cigar is the CIGAR string to parse
#'
#' @return factor of operations
#'
#' @examples
#' cigarStr <-
#'     "46S14M5D14M2D44M3I6M2I14M3I15M1D27M1D11M1I20M1I8M2I2M3I10S"
#' cigarQ2Op(c(20,40,60,80,100), cigarStr)
#'
#' @export
cigarQ2Op <- function(poi, cigar) {
    pmap <- buildPositionMap(cigar)
    return(pmap[match(poi, pmap$qpos),"event"])
}


