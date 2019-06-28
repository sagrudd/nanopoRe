
#' plot a tiled panel of chromosomal depths-of-coverage
#'
#' This method will display depths-of-coverage for the chromosomes contained in the provided parsedBamFile
#' container
#'
#' @param parsedBamFile is the location to the BAM file to parse
#' @param colMax defines the number of columns
#'
#' @examples
#' \dontrun{
#' plotDepthOfCoverageMegablock(parsedBamFile)
#' }
#'
#' @export
plotDepthOfCoverageMegablock <- function(parsedBamFile, colMax=4) {

  suppressWarnings(posMatrix <- matrix(gtools:::mixedsort(unique(parsedBamFile$chrId)), ncol=colMax, byrow=TRUE))
  # data may be recycled ... remove duplicate values ...
  posMatrix[which(duplicated(posMatrix[seq(nrow(posMatrix) * ncol(posMatrix))]))]<-NA
  megaCov.df <- cbind(parsedBamFile,
                      as.data.frame(matrix(unlist(lapply(parsedBamFile$chrId, function(x) { which(posMatrix==x, arr.ind=TRUE) })), ncol=2, byrow=TRUE, dimnames=list(c(), c("row", "col"))))
  )

  meanCov <- mean(parsedBamFile$meanCov, na.rm=TRUE)

  megaCov.df$startpos <- (as.numeric(megaCov.df$startpos - 1))/1000000
  # set ylim to mean + 2stdev?
  ylimit <- mean(megaCov.df$meanCov, na.rm=TRUE) + 2*sd(megaCov.df$meanCov, na.rm=TRUE)

  plotLegend <- paste("chr",gtools:::mixedsort(unique(parsedBamFile$chrId)))
  plotCols <- ceiling(length(plotLegend) / colMax)
  legendDF <- data.frame(x=Inf, y=Inf,
                         lab=plotLegend,
                         row=unlist(lapply(1:plotCols, rep, times=colMax))[1:length(plotLegend)],
                         col=rep(seq(1, colMax), length.out=length(plotLegend)))

  suppressWarnings(megadepthplot <- ggplot(megaCov.df, aes(startpos, meanCov)) +
                     geom_hline(yintercept=meanCov, size=0.3, colour="red") +
                     geom_line(colour=brewer.pal(6, "Paired")[2]) +
                     facet_grid(rows = vars(row), cols=vars(col), shrink=TRUE) +
                     coord_cartesian(ylim=c(0, ylimit))  +
                     scale_x_continuous(labels = comma) +
                     theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text.y = element_blank(), strip.text.x = element_blank()) +
                     xlab("Position along chromosome (Mb)") +
                     ylab("Depth of Coverage (X)") +
                     labs(title="Plot showing depth of coverage vs position for chromosomes mapped") +
                     geom_text(aes(x,y,label=lab), data=legendDF, vjust=1, hjust=1, size=3.5) + theme(plot.title = element_text(size=11)))

  return(megadepthplot)
}

