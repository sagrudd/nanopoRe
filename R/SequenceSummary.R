
#' load a sequencing_summary.txt file into memory
#'
#' importSequencingSummary loads a sequencing_summary.txt file into
#' memory and performs basic sanity checking to ensure that the file
#' is cleaned of potential duplicate headers from aggregation
#'
#' @param seqsum is a path to a file
#' @return data.frame of observations from the sequencing_summary.txt file provided
#'
#' @examples
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#'
#' @export
importSequencingSummary <- function(seqsum) {
  # downsample performed with
  #     cat lambda_sequencing_summary.txt | awk ' BEGIN {srand()} {print rand() " " $0}'  | sort | head -5 | sed 's/[^ ]* //'
  seqsumdata <- data.table::fread(seqsum, stringsAsFactors=FALSE)

  # remove the redundant headers from merged files
  if (length(which(seqsumdata[,1]=="filename")) > 0) {
    seqsumdata <- seqsumdata[-which(seqsumdata[,1]=="filename"),]
  }

  # coerce the columns used in analytics into more appropriate data-types
  seqsumdata$channel<-as.numeric(seqsumdata$channel)
  seqsumdata$start_time<-as.numeric(seqsumdata$start_time)
  seqsumdata$duration<-as.numeric(seqsumdata$duration)
  seqsumdata$num_events<-as.numeric(seqsumdata$num_events)
  seqsumdata$sequence_length_template<-as.numeric(seqsumdata$sequence_length_template)
  seqsumdata$mean_qscore_template<-as.numeric(seqsumdata$mean_qscore_template)

  # passes_filtering is a useful flag; but there are examples of sequencing_summary.txt where this
  # is not present - https://github.com/a-slide/pycoQC/blob/master/pycoQC/data/sequencing_summary_1D_DNA_Albacore_1.2.1.txt
  if (! "passes_filtering" %in% colnames(seqsumdata)) {
    # set all of the reads to pass? apply a cutoff?
    seqsumdata$passes_filtering <- TRUE
  } else {
    seqsumdata$passes_filtering <- as.logical(seqsumdata$passes_filtering)
  }

  setCachedObject("seqsumdata", seqsumdata)

  return(seqsumdata)
}





#' prepare a gauge plot of sequencing_summary reads passing QC
#'
#' plots the basic eye-candy gauge plot of reads passing QC threshold
#'
#' @importFrom dplyr mutate
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @return ggplot2 gauge plot
#'
#' @examples
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- sequencingSummaryPassGauge(seqsum)
#'
#' @export
sequencingSummaryPassGauge <- function(seqsum=NA) {

  if (!is.data.frame(seqsum) && is.na(seqsum)) {
    oname <- "seqsumdata"
    if (hasCachedObject(oname)) {
      seqsum <- getCachedObject(oname)
    }
  }

  df <- data.frame(matrix(nrow=1, ncol = 3))

  names(df) <- c("variable", "percentage","label")
  df$variable <- c("pass")
  df$percentage <- c(round(length(which(seqsum$passes_filtering==TRUE)) / nrow(seqsum), 3))

  df <- df %>% dplyr::mutate(group=ifelse(percentage <0.6, "red",
                                   ifelse(percentage>=0.6 & percentage<0.8, "orange","green")),
                      label=paste0(df$percentage*100, "%"))

  title="Percentage of reads\npassing QC filter"

  gaugePlot <- ggplot(df, aes(fill = group, ymax = percentage, ymin = 0, xmax = 2, xmin = 1)) +
    geom_rect(aes(ymax=1, ymin=0, xmax=2, xmin=1), fill ="#ece8bd") +
    geom_rect() +
    coord_polar(theta = "y",start=-pi/2) + xlim(c(0, 2)) + ylim(c(0,2)) +
    guides(fill=FALSE) +
    guides(colour=FALSE) +
    theme_void() +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    geom_text(aes(x = 0, y = 0, label = label), size=13) +
    geom_text(aes(x=1.5, y=1.5, label=title), size=11) +
    scale_fill_manual(values = c("red"="#C9146C", "orange"="#DA9112", "green"="#129188")) +
    scale_colour_manual(values = c("red"="#C9146C", "orange"="#DA9112", "green"="#129188"))

  return(gaugePlot)
}






#' prepare a gauge plot of sequencing_summary reads passing QC
#'
#' plots the basic eye-candy gauge plot of reads passing QC threshold
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param platform is the nanopore platform [MinION/Flongle/PromethION]
#' @return ggplot2 gauge plot
#'
#' @examples
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- sequencingSummaryChannelActivity(seqsum)
#'
#' @export
sequencingSummaryChannelActivity <- function(seqsum=NA, platform=NA) {

  if (!is.data.frame(seqsum) && is.na(seqsum)) {
    oname <- "seqsumdata"
    if (hasCachedObject(oname)) {
      seqsum <- getCachedObject(oname)
    }
  }

  if (is.na(platform)) {
    platform <- sequencingSummaryGetPlatform(seqsum)
  }

  channelMap <- sequencingSummaryGetChannelMap(platform)


}


#' identify the most likely sequencing platform used to create the summary data
#'
#' when provided with the seqsum object from a Sequencing_summary.txt file identify the
#' most likely sequencing platform used
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @return platform is the nanopore platform [MinION/Flongle/PromethION]
#'
#' @examples
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' platform <- sequencingSummaryGetPlatform(seqsum)
#'
#' @export
sequencingSummaryGetPlatform <- function(seqsum=NA) {

  if (!is.data.frame(seqsum) && is.na(seqsum)) {
    oname <- "seqsumdata"
    if (hasCachedObject(oname)) {
      seqsum <- getCachedObject(oname)
    }
  }

  platform <- "MinION"

  if (max(seqsum$channel) < 130) {
    # this is likely to be a Flongle ...
    platform <- "Flongle"
  }

  if (max(seqsum$channel) > 1000) {
    # this is likely to be a PromethION
    platform <- "PromethION"
  }

  return(platform)
}






sequencingSummaryGetChannelMap <- function(platform) {
  if (platform=="MinION") {
    return(getMinIONChannelMap())
  } else if (plaform=="Flongle") {
    return(getFlongleChannelMap())
  } else if (platform=="PromethION") {
    return(getPromethIONChannelMap())
  }
  return(NULL)
}



#' produce the channelMap for a MinION flowcell for spatial plots
#'
#' prepares a data.frame of channelIds and their X, Y coordinates for a MinION flowcell
#'
#' @return data.frame with channel, row and col columns
#'
#' @examples
#' getMinIONChannelMap()
#'
#' @export
getMinIONChannelMap <- function() {
  # build the map for R9.4.1 flowcell, as a long-form dataframe
  blockCalc <- function(i) {
    m <- matrix(seq(i,i+63,by=1), ncol=8, byrow=TRUE)
    cbind(m[seq(5,8,by=1), ], m[seq(4), rev(seq(8))])
  }
  layout <- do.call(rbind, lapply(c(1, 449, 385, 321, 257, 193, 129, 65), blockCalc))
  channelMap <- as.data.frame(cbind(channel=as.vector(t(layout)), which(layout==as.vector(layout), arr.ind=TRUE)))
  return(channelMap)
}


#' produce the channelMap for a flongle flowcell for spatial plots
#'
#' prepares a data.frame of channelIds and their X, Y coordinates for a flongle flowcell
#'
#' @return data.frame with channel, row and col columns
#'
#' @examples
#' getFlongleChannelMap()
#'
#' @export
getFlongleChannelMap <- function() {
  layout <- matrix(c(seq(1,12), 0, seq(13,24), 0, seq(25,114), 0, seq(115, 126), 0), ncol=13, byrow=TRUE)
  layout <- layout[rev(seq(10)),]
  channelMap <- as.data.frame(cbind(channel=as.vector(t(layout)), which(layout==as.vector(layout), arr.ind=TRUE)))
  return(channelMap)
}


getPromethIONChannelMap <- function() {

}

