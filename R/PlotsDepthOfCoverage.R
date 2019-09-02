
#' plot a tiled panel of chromosomal depths-of-coverage
#'
#' This method will display depths-of-coverage for the chromosomes contained in the provided parsedBamFile
#' container
#'
#' @import ggplot2
#' @importFrom stats sd
#' @import RColorBrewer
#' @return ggplot2 image
#' @param coverageData data.frame of coverage data
#' @param colMax defines the number of columns
#'
#' @examples
#' \dontrun{
#' coverageData <- bamSummaryToCoverage(bamFile, tilewidth=250000)
#' plotDepthOfCoverageMegablock(parsedBamFile)
#' }
#'
#' @export
plotDepthOfCoverageMegablock <- function(coverageData, colMax = 4) {
    
    coverageDF <- as.data.frame(coverageData, stringsAsFactors = FALSE)
    coverageDF$seqnames <- as.character(coverageDF$seqnames)
    
    suppressWarnings(posMatrix <- matrix(gtools::mixedsort(unique(coverageDF$seqnames)), ncol = colMax, 
        byrow = TRUE))
    # data may be recycled ... remove duplicate values ...
    posMatrix[which(duplicated(posMatrix[seq(nrow(posMatrix) * ncol(posMatrix))]))] <- NA
    coverageDF <- cbind(coverageDF, as.data.frame(matrix(unlist(lapply(coverageDF$seqnames, function(x) {
        which(posMatrix == x, arr.ind = TRUE)
    })), ncol = 2, byrow = TRUE, dimnames = list(c(), c("row", "col")))))
    
    meanCov <- mean(coverageDF$binned_cov, na.rm = TRUE)
    
    coverageDF$start <- (as.numeric(coverageDF$start - 1))/1e+06
    # set ylim to mean + 2stdev?
    ylimit <- mean(coverageDF$binned_cov, na.rm = TRUE) + 2 * stats::sd(coverageDF$binned_cov, na.rm = TRUE)
    
    plotLegend <- paste("chr", gtools::mixedsort(unique(coverageDF$seqnames)))
    plotCols <- ceiling(length(plotLegend)/colMax)
    legendDF <- data.frame(x = Inf, y = Inf, lab = plotLegend, row = unlist(lapply(seq_len(plotCols), 
        rep, times = colMax))[seq_along(plotLegend)], col = rep(seq(1, colMax), length.out = length(plotLegend)))
    
    suppressWarnings(megadepthplot <- ggplot(coverageDF, aes_string("start", "binned_cov")) + geom_hline(yintercept = meanCov, 
        size = 0.3, colour = "red") + geom_line(colour = brewer.pal(6, "Paired")[2]) + facet_grid(rows = vars(row), 
        cols = vars(col), shrink = TRUE) + coord_cartesian(ylim = c(0, ylimit)) + scale_x_continuous(labels = comma) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text.y = element_blank(), strip.text.x = element_blank()) + 
        xlab("Position along chromosome (Mb)") + ylab("Depth of Coverage (X)") + labs(title = "Plot showing depth of coverage vs position for chromosomes mapped") + 
        geom_text(aes_string("x", "y", label = "lab"), data = legendDF, vjust = 1, hjust = 1, size = 3.5) + 
        theme(plot.title = element_text(size = 11)))
    
    return(megadepthplot)
}




#' plot a histogram of whole genome depth-of-coverage
#'
#' This method will plot a histogram of whole genome depth-of-coverage
#'
#' @importFrom utils head
#' @param bamFile is the path to the bamFile to plot
#' @param segments number of segments to split data into
#' @return ggplot2 plot
#'
#' @examples
#' \dontrun{
#' plotDepthOfCoverageMegablock(bamFile)
#' }
#'
#' @export
plotOverallCovHistogram <- function(bamFile, segments = 50) {
    
    coverage <- bamSummaryToCoverage(bamFile)
    
    meanCov <- mean(coverage$binned_cov, na.rm = TRUE)
    
    coverageMatrix <- data.frame(table(coverage$binned_cov))
    coverageMatrix$bases <- as.numeric(coverageMatrix$Freq * mean(width(coverage)))
    coverageMatrix[, 1] <- as.numeric(levels(coverageMatrix[, 1]))[coverageMatrix[, 1]]
    # mask the regions with unmapped reads
    coverageMatrix <- coverageMatrix[-which(coverageMatrix[, 1] == 0), ]
    
    breaks <- seq(0, 2 * meanCov, length.out = segments)
    binAssignments <- cut(coverageMatrix$Var1, breaks, include.lowest = TRUE, right = FALSE)
    
    scrapeBinnedBases <- function(level) {
        sum(coverageMatrix[which(binAssignments == level), "bases"])
    }
    binnedData <- unlist(lapply(levels(binAssignments), scrapeBinnedBases))
    binnedDf <- data.frame(coverage = utils::head(breaks, -1), bases = binnedData)
    
    
    plot <- ggplot(binnedDf, aes_string(x = "coverage", y = "bases")) + geom_vline(xintercept = meanCov, 
        size = 0.3, colour = "red") + geom_bar(stat = "identity", fill = brewer.pal(6, "Paired")[2]) + 
        scale_y_continuous(labels = comma) + labs(title = "Histogram showing distribution of depth-of-coverage across genome") + 
        xlab("Coverage (x)") + ylab("Number of bases") + theme(plot.title = element_text(size = 11))
    
    return(plot)
}

