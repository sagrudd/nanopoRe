

#' plot a histogram of mapping accuracies
#'
#' This method will plot a histogram of mapping accuracies
#'
#' @usage plotAlignmentAccuracy(bamFile, flag = "Primary", lower = 0.7,
#'     upper = 1, segments = 50)
#' @param bamFile is the path to the bamFile to plot
#' @param flag is the flagged sequence type to plot [Primary]
#' @param lower lower limit for plotting values
#' @param upper upper limit for plotting values
#' @param segments number of segments to split data into
#' @return a ggplot2 object
#'
#' @examples
#' \dontrun{
#' plotAlignmentAccuracy(bamFile)
#' }
#'
#' @export
plotAlignmentAccuracy <- function(bamFile, flag = "Primary", lower = 0.7, upper = 1, segments = 50) {

    alignmentData <- bamSummarise(bamFile)
    # slice data with dplyr
    filteredData <- alignmentData %>% filter(.data$readFlag == flag & .data$accuracy >= lower & .data$accuracy <= upper)

    meanAccuracy <- mean(filteredData$accuracy, na.rm = TRUE) * 100

    accuracyMatrix <- data.frame(table(round(filteredData$accuracy * 100, digits = 1)))
    accuracyMatrix[, 1] <- as.numeric(levels(accuracyMatrix[, 1]))[accuracyMatrix[, 1]]

    breaks <- seq(lower * 100, (upper * 100) + 1, length.out = segments)
    binAssignments <- cut(accuracyMatrix$Var1, breaks, include.lowest = TRUE, right = FALSE)

    scrapeBinnedBases <- function(level) {
        sum(accuracyMatrix[which(binAssignments == level), "Freq"])
    }
    binnedData <- unlist(lapply(levels(binAssignments), scrapeBinnedBases))
    binnedDf <- data.frame(accuracy = utils::head(breaks, -1), reads = binnedData)

    plot <- ggplot(binnedDf, aes_string(x = "accuracy", y = "reads")) + geom_vline(xintercept = meanAccuracy, size = 0.3, colour = "red") +
        geom_bar(stat = "identity", fill = brewer.pal(6, "Paired")[2]) + scale_y_continuous(labels = comma) + labs(title = "Histogram showing distribution of mapping accuracies") +
        xlab("Accuracy (%)") + ylab("Number of reads") + theme(plot.title = element_text(size = 11))

    return(plot)
}



#' plot a histogram of mapping identities
#'
#' This method will plot a histogram of mapping identities
#'
#' @usage plotAlignmentIdentity(bamFile, flag = "Primary", lower = 0.7,
#'     upper = 1, segments = 50)
#' @param bamFile is the path to the bamFile to plot
#' @param flag is the flagged sequence type to plot [Primary]
#' @param lower lower limit for plotting values
#' @param upper upper limit for plotting values
#' @param segments number of segments to split data into
#' @return a ggplot2 object
#'
#' @examples
#' \dontrun{
#' plotAlignmentIdentity(bamFile)
#' }
#'
#' @export
plotAlignmentIdentity <- function(bamFile, flag = "Primary", lower = 0.7, upper = 1, segments = 50) {

    alignmentData <- bamSummarise(bamFile)
    # slice data with dplyr
    filteredData <- alignmentData %>% filter(.data$readFlag == flag & .data$identity >= lower & .data$identity <= upper)

    meanAccuracy <- mean(filteredData$identity, na.rm = TRUE) * 100

    identityMatrix <- data.frame(table(round(filteredData$identity * 100, digits = 1)))
    identityMatrix[, 1] <- as.numeric(levels(identityMatrix[, 1]))[identityMatrix[, 1]]

    breaks <- seq(lower * 100, (upper * 100) + 1, length.out = segments)
    binAssignments <- cut(identityMatrix$Var1, breaks, include.lowest = TRUE, right = FALSE)

    scrapeBinnedBases <- function(level) {
        sum(identityMatrix[which(binAssignments == level), "Freq"])
    }
    binnedData <- unlist(lapply(levels(binAssignments), scrapeBinnedBases))
    binnedDf <- data.frame(identity = utils::head(breaks, -1), reads = binnedData)

    plot <- ggplot(binnedDf, aes_string(x = "identity", y = "reads")) + geom_vline(xintercept = meanAccuracy, size = 0.3, colour = "red") +
        geom_bar(stat = "identity", fill = brewer.pal(6, "Paired")[2]) + scale_y_continuous(labels = comma) + labs(title = "Histogram showing distribution of mapping identities") +
        xlab("Identity (%)") + ylab("Number of reads") + theme(plot.title = element_text(size = 11))

    return(plot)
}



#' plot two ggplot2 figures side by side
#'
#' This method will plot 2 figures in a single plot space
#'
#' @import grid
#' @param plotLeft is the left-hand plot
#' @param plotRight is the right-hand plot
#' @param ldim relative scaling for left plot (2)
#' @param rdim relative scaling for right plot (2)
#' @return a ggplot2 object
#'
#' @examples
#' \dontrun{
#' LeftRightPlot(p1, p2)
#' }
#'
#' @export
LeftRightPlot <- function(plotLeft, plotRight, ldim = 2, rdim = 2) {
    grid.newpage()
    # note that this is not x and y but rows and columns
    pushViewport(viewport(layout = grid.layout(1, (ldim + rdim))))
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    print(plotLeft, vp = vplayout(1, seq(ldim)))
    print(plotRight, vp = vplayout(1, seq(rdim) + ldim))
}


