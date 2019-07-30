
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
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' }
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
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- sequencingSummaryPassGauge(seqsum)
#' }
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






#' prepare a channel activity plot from sequencing_summary reads file
#'
#' plots the basic eye-candy gauge channel activity plot of reads against channel of origin
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param platform is the nanopore platform [MinION/Flongle/PromethION]
#' @return ggplot2 channel activity plot
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- sequencingSummaryChannelActivity(seqsum)
#' }
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


  hm.palette <- colorRampPalette(brewer.pal(9, 'Blues'), space='Lab') #RdPu, Oranges, Greens, YlOrRd, Purples

  channelCounts <- as.data.frame(matrix(rep(0, max(channelMap$channel)), ncol=1))
  channelCountRaw <- as.data.frame(table(unlist(seqsum[, "channel"])), row.names=1)
  channelCounts[row.names(channelCountRaw),] <- channelCountRaw[,1]

  channelMap <- merge(channelMap, channelCounts, by.x="channel", by.y=0)
  colnames(channelMap)[4]<-"count"
  channelMapMatrix <- reshape2::acast(channelMap, col ~ row, value.var = "count")

  theme_update(plot.title = element_text(hjust = 0.5))

  activityPlot <- ggplot(channelMap, aes(x = row, y = col, fill = count)) +
    geom_tile() +
    geom_text(data=channelMap,aes(x=row, y=col,label=count,color=count),show.legend = F, size=2.5) +
    scale_x_discrete(breaks=NULL) +
    scale_y_discrete(breaks=NULL) +
    coord_equal() +
    scale_fill_gradientn(colours = hm.palette(100)) +
    scale_color_gradient2(low = hm.palette(100), high = hm.palette(1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="Channel activity plot showing number of reads per flowcell channel") +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="bottom",
          legend.key.width=unit(5.6,"cm"))

  return(activityPlot)
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
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' platform <- sequencingSummaryGetPlatform(seqsum)
#' }
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







# https://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x
roundUpNice <- function(x, nice=seq(from=1, to=10, by=0.25)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}


getBinAssignments <- function(seqsum, breaks) {
  binAssignments <- cut(seqsum$sequence_length_template, breaks, include.lowest=TRUE, right=FALSE)
  return(binAssignments)
}

#' @importFrom stats quantile
getBinBreaks <- function(seqsum) {
  # pick a friendly upper limit to render sequence lengths into a histogram
  # here we're aiming for a robustly rounded up 97.5 quantile of the data (skip a few outliers ...)
  upperLimit <- roundUpNice(as.numeric(quantile(x=seqsum$sequence_length_template, probs=c(0.975))))
  # an ideal histogram will have 40 or so bins
  histogramBinCount <- 40
  breakVal = roundUpNice(upperLimit / histogramBinCount)
  breaks <- seq(0, to=upperLimit, by=breakVal)
  return(breaks)
}


#' plot a weighted histogram of sequence read lengths
#'
#' plots the histogram of read lengths, weighted, and shaded by pass/fail status
#'
#' @importFrom dplyr last
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @return ggplot2 showing weighted read length distribution
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- sequencingSummaryWeightedReadLength(seqsum)
#' }
#'
#' @export
sequencingSummaryWeightedReadLength <- function(seqsum) {

  breaks <- getBinBreaks(seqsum)
  breakVal <- breaks[2] # assuming that the range is 0 based
  upperLimit <- dplyr::last(breaks)
  binAssignments <- getBinAssignments(seqsum, breaks)

  passedSeqs <- seqsum[which(seqsum$passes_filtering), ]
  N50 <- ncalc(passedSeqs$sequence_length_template, 0.5)
  passedMeanLength = round(mean(passedSeqs$sequence_length_template), digits = 0)

  scrapeBinnedBases <- function(level, qcpass, binAssignments, seqsum) {
    sum(subset(seqsum[which(binAssignments == level), ], passes_filtering==qcpass)$sequence_length_template)
  }

  passedBinnedBases <- unlist(lapply(levels(binAssignments), scrapeBinnedBases, qcpass=TRUE, binAssignments=binAssignments, seqsum=seqsum))
  failedBinnedBases <- unlist(lapply(levels(binAssignments), scrapeBinnedBases, qcpass=FALSE, binAssignments=binAssignments, seqsum=seqsum))

  binnedBaseDist <- data.frame(length=head(breaks, -1), pass=passedBinnedBases, fail=failedBinnedBases)
  binnedBaseMelt <- reshape2::melt(binnedBaseDist, id.vars=c("length"))

  weightedReadLengths <- ggplot(binnedBaseMelt, aes_string(x="length", fill="variable", y="value")) +
    geom_bar(stat="identity") +
    xlab("Read length\n") + ylab("Number of bases sequenced\n") +
    scale_fill_manual("QC", values=c("fail"=brewer.pal(6, "Paired")[1], "pass"=brewer.pal(6, "Paired")[2])) +
    scale_x_continuous(limits=c(-breakVal,upperLimit), breaks=pretty(passedSeqs$sequence_length_template,n=40)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="Histogram showing the number of sequenced bases against sequence length", fill="QV filter")+
    geom_vline(xintercept = N50, size = 1) +
    annotate("text", x=N50, y=max(passedBinnedBases + failedBinnedBases), label = " N50", hjust=0, colour="SteelBlue") +
    geom_vline(xintercept = passedMeanLength, size = 1) +
    annotate("text", x=passedMeanLength, y=max(passedBinnedBases + failedBinnedBases), label = " Mean", hjust=0, colour="SteelBlue")

  return(weightedReadLengths)
}







#' plot a histogram of sequence read lengths
#'
#' plots the histogram of read lengths shaded by pass/fail status
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @return ggplot2 showing read length distribution
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- sequencingSummaryReadLengthHistogram(seqsum)
#' }
#'
#' @export
sequencingSummaryReadLengthHistogram <- function(seqsum) {

  breaks <- getBinBreaks(seqsum)
  breakVal <- breaks[2] # assuming that the range is 0 based
  upperLimit <- dplyr::last(breaks)
  binAssignments <- getBinAssignments(seqsum, breaks)

  passedSeqs <- seqsum[which(seqsum$passes_filtering), ]
  N50 <- ncalc(passedSeqs$sequence_length_template, 0.5)
  passedMeanLength = round(mean(passedSeqs$sequence_length_template), digits = 0)

  scrapeBinnedReads <- function(level, qcpass) {
    length(subset(seqsum[which(binAssignments == level), ], passes_filtering==qcpass)$sequence_length_template)
  }

  passedBinnedReads <- unlist(lapply(levels(binAssignments), scrapeBinnedReads, qcpass=TRUE))
  failedBinnedReads <- unlist(lapply(levels(binAssignments), scrapeBinnedReads, qcpass=FALSE))

  binnedReadDist <- data.frame(length=head(breaks, -1), pass=passedBinnedReads, fail=failedBinnedReads)
  binnedReadMelt <- reshape2::melt(binnedReadDist, id.vars=c("length"))

  lengthHistogram <- ggplot(binnedReadMelt, aes(x=length, fill=variable, y=value)) +
    geom_bar(stat="identity") +
    xlab("Read length\n") + ylab("Number of reads\n") +
    scale_fill_manual("QC", values=c("fail"=brewer.pal(6, "Paired")[1], "pass"=brewer.pal(6, "Paired")[2])) +
    scale_x_continuous(limits=c(-breakVal,upperLimit), breaks=pretty(passedSeqs$sequence_length_template,n=40)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="Histogram showing distribution of read lengths across quality passing sequences", fill="QV filter")+
    geom_vline(xintercept = N50, size = 1) +
    annotate("text", x=N50, y=max(passedBinnedReads + failedBinnedReads), label = " N50", hjust=0, colour="SteelBlue") +
    geom_vline(xintercept = passedMeanLength, size = 1) +
    annotate("text", x=passedMeanLength, y=max(passedBinnedReads + failedBinnedReads), label = " Mean", hjust=0, colour="SteelBlue")

  return(lengthHistogram)
}




#' plot a histogram of sequence quality scores
#'
#' plots the histogram of read mean quality scores shaded by pass/fail status
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @return ggplot2 showing read quality distribution
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- sequencingSummaryReadQualityHistogram(seqsum)
#' }
#'
#' @export
sequencingSummaryReadQualityHistogram <- function(seqsum) {
  qdist <- ggplot(seqsum, aes(x=mean_qscore_template, fill=passes_filtering)) +
    geom_histogram(breaks=seq(from=0, to=15, by=0.1)) +
    scale_fill_manual(name="QC", values=c("TRUE"=brewer.pal(6, "Paired")[2], "FALSE"=brewer.pal(6, "Paired")[1]), labels=c( "pass", "fail"), breaks=c("TRUE", "FALSE")) +
    labs(title="Plot showing distribution of quality scores across all reads") +
    xlab("Mean Q score of read") +
    ylab("Number of reads")

  return(qdist)
}






#' plot a density map of sequence lengths and quality scores
#'
#' plots the density plot of read length against read mean quality scores
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param binFilter is the minimum number of reads to include a cell in plot (removes speckle)
#' @param qcThreshold is the QC threshold used
#' @return ggplot2 showing densities of read length and quality distribution
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- sequencingSummaryReadLengthQualityDensity(seqsum)
#' }
#'
#' @export
sequencingSummaryReadLengthQualityDensity <- function(seqsum, binFilter=5, qcThreshold=7) {
  # prepare the density plot, but do not render
  lq_dens <- ggplot(seqsum, aes(log10(sequence_length_template), mean_qscore_template)) + geom_bin2d(bins=100)
  # extract the density map from the plot
  lq_dens_counts <- ggplot_build(lq_dens)$data[[1]]
  if (binFilter > 0) {
    # remove the bins from the density map that do not contain sequence count above threshold
    lq_dens_counts <- lq_dens_counts[-which(lq_dens_counts$count <= binFilter),]
  }
  # directly plot this modified density map (stat=="identity")
  qldensityplot <- ggplot(lq_dens_counts) +
    geom_bin2d(aes(x,y,fill=count), stat="identity") +
    scale_fill_distiller(palette="Blues", trans="reverse") +
    geom_hline(yintercept = qcThreshold, size = 1) +
    xlab("log10(read length)") +
    ylab("read mean quality") +
    scale_x_continuous(breaks = c(1,2,3,4,5), labels = c("10", "100", "1000", "10,000", "100,000")) +
    annotation_logticks(base = 10, sides = "b", scaled = TRUE) +
    labs(title="Contour Plot showing distribution of quality scores against log10 read lengths (all reads)")
  return(qldensityplot)
}



getTemporalDataset <- function(seqsum, sampleIntervalMinutes, breaks, binass) {
  mergeItPerHour <- function(interval, binnedAssignments, filter) {
    totalbases = 0
    if (length(which(binnedAssignments==interval))>0) {
      subset <- seqsum[which(binnedAssignments==interval), ]
      if (length(which(subset$passes_filtering == filter)) > 0) {
        totalbases = sum(subset[which(subset$passes_filtering == filter), "sequence_length_template"])
      }
    }
    # need to scale what is being returned - totalbases value is total bases within an interval (sampleIntervalMinutes)
    return(totalbases / 1e9 / sampleIntervalMinutes * 60)
  }

  binnedTemporalDataPerHour <- data.frame(
    cbind(
      time=breaks,
      pass=unlist(lapply(seq(breaks), mergeItPerHour, binnedAssignments=binass,filter=TRUE)),
      fail=unlist(lapply(seq(breaks), mergeItPerHour, binnedAssignments=binass, filter=FALSE))
    )
  )

  binnedTemporalDataPerHour$time <- binnedTemporalDataPerHour$time / 60 / 60
  return(binnedTemporalDataPerHour)
}



#' plot a sequence throughput against time for specified sequencing_summary run
#'
#' plots a ggplot2 graph of performance against time for run and separates passed and failed sequence reads
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param scaling scale factor for the data
#' @param sampleHours is the number of hours to plot data for (default is 48)
#' @param sampleIntervalMinutes is the resolution to plot data at
#' @return ggplot2 showing temporal performance
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- SequencingSummaryTemporalThroughput(seqsum)
#' }
#'
#' @export
SequencingSummaryTemporalThroughput <- function(seqsum, scaling=1, sampleHours = 48, sampleIntervalMinutes = 60) {

  seqsum$start_time <- seqsum$start_time - min(seqsum$start_time)
  seqsum$start_time <- seqsum$start_time / scaling

  breaks = seq(0, sampleHours*60*60, by=60*sampleIntervalMinutes)
  binass <- findInterval(seqsum$start_time, breaks)

  binnedTemporalDataPerHour <- getTemporalDataset(seqsum, sampleIntervalMinutes, breaks, binass)

  plot <- ggplot(binnedTemporalDataPerHour, aes(time)) +
    geom_line(aes(y = fail, colour = "fail"), size=1) +
    geom_line(aes(y = pass, colour = "pass"), size=1) +
    scale_color_manual(name="QV", values=c("fail"=brewer.pal(6, "Paired")[1], "pass"=brewer.pal(6, "Paired")[2])) +
    xlab("Time (hours)") +
    ylab("Gigabases sequenced per hour") +
    labs(title="Plot showing sequence throughput against time")

  return(plot)
}


#' plot cumulative volumes of sequence bases
#'
#' plots a ggplot2 graph of accumulated sequenced bases against time for run and separates bases called
#' from passed and failed sequence reads
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param scaling scale factor for the data
#' @param sampleHours is the number of hours to plot data for (default is 48)
#' @param sampleIntervalMinutes is the resolution to plot data at
#' @return ggplot2 showing temporal performance
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- SequencingSummaryCumulativeBases(seqsum)
#' }
#'
#' @export
SequencingSummaryCumulativeBases <- function(seqsum, scaling=1, sampleHours = 48, sampleIntervalMinutes = 60) {

  seqsum$start_time <- seqsum$start_time - min(seqsum$start_time)
  seqsum$start_time <- seqsum$start_time / scaling

  breaks = seq(0, sampleHours*60*60, by=60*sampleIntervalMinutes)
  binass <- findInterval(seqsum$start_time, breaks)

  binnedTemporalDataPerHour <- getTemporalDataset(seqsum, sampleIntervalMinutes, breaks, binass)

  # binnedTemporalDataPerHour is scaled to Gbp per hour - rescale to raw for cumulative plotting
  binnedTemporalDataPerHour$pass <- binnedTemporalDataPerHour$pass / 60 * sampleIntervalMinutes
  binnedTemporalDataPerHour$fail <- binnedTemporalDataPerHour$fail / 60 * sampleIntervalMinutes

  base50 <- SequencingSummaryBase50(seqsum, b=0.5)
  base90 <- SequencingSummaryBase50(seqsum, b=0.9)
  T50 <- SequencingSummaryT50(seqsum, t=0.5, scaling=scaling, sampleHours=sampleHours, sampleIntervalMinutes=sampleIntervalMinutes)
  T90 <- SequencingSummaryT50(seqsum, t=0.9, scaling=scaling, sampleHours=sampleHours, sampleIntervalMinutes=sampleIntervalMinutes)


  cumulativePlot <- ggplot(binnedTemporalDataPerHour, aes(time)) +
    geom_line(aes(y = cumsum(fail), colour = "fail"), size=1) +
    geom_line(aes(y = cumsum(pass), colour = "pass"), size=1) +
    scale_color_manual(name="QV", values=c("fail"=brewer.pal(6, "Paired")[1], "pass"=brewer.pal(6, "Paired")[2])) +
    geom_segment(x=T50$minimum, y=0, xend=T50$minimum, yend=base50, colour="darkgray", size=1) +
    geom_segment(x=0, y=base50, xend=T50$minimum, yend=base50, colour="darkgray", size=1) +
    annotate("text", x=T50$minimum, y=base50, label=" T50", vjust=1, hjust=0, colour="SteelBlue") +
    geom_segment(x=T90$minimum, y=0, xend=T90$minimum, yend=base90, colour="darkgray", size=1) +
    geom_segment(x=0, y=base90, xend=T90$minimum, yend=base90, colour="darkgray", size=1) +
    annotate("text", x=T90$minimum, y=base90, label=" T90", vjust=1, hjust=0, colour="SteelBlue") +
    xlab("Time (hours)") +
    ylab("Number of bases sequenced (Gigabases)") +
    labs(title="Plot showing cumulative bases sequenced against time")

  return(cumulativePlot)

}




#' calculates the fractional number of bases according to supplied b parameter
#'
#' an accessory method for various logicals; simple fractional base calculator
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param b is a fractional point through run against which time will be calculated
#' @return a numeric value expressed in gigabases
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' Base50 <- SequencingSummaryBase50(seqsum)
#' }
#'
#' @export
SequencingSummaryBase50 <- function(seqsum, b=0.5) {
  passedSeqs <- seqsum[which(seqsum$passes_filtering), ]
  base50 <- sum(passedSeqs$sequence_length_template) / 1e9 * b
  return(base50)
}




#' calculates the timepoint within a sequencing run where 50percent of the data is produced
#'
#' an accessory method for identifying a timepoint where a given amount of data has been produced
#'
#' @importFrom stats approxfun
#' @importFrom stats optimize
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param t is a fractional point through run against which time will be calculated
#' @param scaling factor
#' @param sampleHours is the number of hours to consider
#' @param sampleIntervalMinutes is the resolution of the plot in minutes
#' @return a numeric value expressed in hours
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' T50 <- SequencingSummaryT50(seqsum)
#' }
#'
#' @export
SequencingSummaryT50 <- function(seqsum, t=0.5, scaling=1, sampleHours = 48, sampleIntervalMinutes = 60) {

  seqsum$start_time <- seqsum$start_time - min(seqsum$start_time)
  seqsum$start_time <- seqsum$start_time / scaling

  breaks = seq(0, sampleHours*60*60, by=60*sampleIntervalMinutes)
  binass <- findInterval(seqsum$start_time, breaks)

  binnedTemporalDataPerHour <- getTemporalDataset(seqsum, sampleIntervalMinutes, breaks, binass)

  # binnedTemporalDataPerHour is scaled to Gbp per hour - rescale to raw for cumulative plotting
  binnedTemporalDataPerHour$pass <- binnedTemporalDataPerHour$pass / 60 * sampleIntervalMinutes

  # https://stackoverflow.com/questions/31404679/can-ggplot2-find-the-intersections-or-is-there-any-other-neat-way
  acquireTimePoints <- which(binnedTemporalDataPerHour$pass > 0)
  targetInterpolate <- approxfun(x=binnedTemporalDataPerHour[acquireTimePoints, "time"], y=cumsum(binnedTemporalDataPerHour[acquireTimePoints, "pass"]))

  base50 <- SequencingSummaryBase50(seqsum, b=t)
  T50 <- optimize(function(t0) abs(targetInterpolate(t0) - base50),
                  interval = range(binnedTemporalDataPerHour[acquireTimePoints, "time"]))

  return(T50)
}




#' plot cumulative volumes of sequence reads
#'
#' plots a ggplot2 graph of accumulated sequence reads against time for run and separates passed and
#' failed sequence reads
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param scaling scale factor for the data
#' @param sampleHours is the number of hours to plot data for (default is 48)
#' @param sampleIntervalMinutes is the resolution to plot data at
#' @return ggplot2 showing temporal performance
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- SequencingSummaryCumulativeReads(seqsum)
#' }
#'
#' @export
SequencingSummaryCumulativeReads <- function(seqsum, scaling=1, sampleHours = 48, sampleIntervalMinutes = 60) {

  seqsum$start_time <- seqsum$start_time - min(seqsum$start_time)
  seqsum$start_time <- seqsum$start_time / scaling

  breaks = seq(0, sampleHours*60*60, by=60*sampleIntervalMinutes)
  binass <- findInterval(seqsum$start_time, breaks)

  mergeItReadsPerHour <- function(interval, binnedAssignments,filter) {
    totalreads = 0
    if (length(which(binnedAssignments==interval))>0) {
      subset <- seqsum[which(binnedAssignments==interval), ]
      if (length(which(subset$passes_filtering == filter)) > 0) {
        totalreads = nrow(subset[which(subset$passes_filtering == filter),])
      }
    }
    # scale results to mean millions of reads per hour
    return(totalreads/ 1e6 / sampleIntervalMinutes * 60)
  }

  binnedTemporalDataReadsPerHour <- data.frame(
    cbind(time=breaks,
          pass=unlist(lapply(seq(breaks), mergeItReadsPerHour, binnedAssignments=binass, filter=TRUE)),
          fail=unlist(lapply(seq(breaks), mergeItReadsPerHour, binnedAssignments=binass, filter=FALSE))
    )
  )

  binnedTemporalDataReadsPerHour$time <- binnedTemporalDataReadsPerHour$time / 60 / 60
  # binnedTemporalDataReadsPerHour is scaled to Gbp per hour - rescale to raw for cumulative plotting
  binnedTemporalDataReadsPerHour$pass <- binnedTemporalDataReadsPerHour$pass / 60 * sampleIntervalMinutes
  binnedTemporalDataReadsPerHour$fail <- binnedTemporalDataReadsPerHour$fail / 60 * sampleIntervalMinutes

  cumulativePlot <- ggplot(binnedTemporalDataReadsPerHour, aes(time)) +
    geom_line(aes(y = cumsum(fail), colour = "fail"), size=1) +
    geom_line(aes(y = cumsum(pass), colour = "pass"), size=1) +
    scale_color_manual(name="QV", values=c("fail"=brewer.pal(6, "Paired")[1], "pass"=brewer.pal(6, "Paired")[2])) +
    xlab("Time (hours)") +
    ylab("Number of reads sequenced (Millions)") +
    labs(title="Plot showing cumulative reads sequenced against time")

  return(cumulativePlot)
}





#' plot speed of sequencing against time (bases per second distribution)
#'
#' plots a ggplot2 box-and-whisker plot for the distribution of sequencing speeds against time
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param scaling scale factor for the data
#' @param sampleHours is the number of hours to plot data for (default is 48)
#' @param sampleIntervalMinutes is the resolution to plot data at
#' @return ggplot2 showing temporal performance
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- SequencingSummarySpeedPlot(seqsum)
#' }
#'
#' @export
SequencingSummarySpeedPlot <- function(seqsum, scaling=1, sampleHours = 48, sampleIntervalMinutes = 60) {

  seqsum$start_time <- seqsum$start_time - min(seqsum$start_time)
  seqsum$start_time <- seqsum$start_time / scaling

  breaks = seq(0, sampleHours*60*60, by=60*sampleIntervalMinutes)
  binass <- findInterval(seqsum$start_time, breaks)

  speedTime <- data.frame(segment=binass, rate=seqsum$sequence_length_template / (seqsum$duration/scaling))

  speedplot <- ggplot(speedTime, aes(x=segment, y=rate, group=segment)) +
    geom_boxplot(fill="steelblue", outlier.shape=NA) +
    scale_x_continuous(name="Time (hours)") +
    ylab("Sequencing rate (bases per second)") +
    labs(title="boxplot showing distribution of translocation speed against time")

  return(speedplot)
}




#' plot number of observed channels actively producing data against time
#'
#' plots a ggplot2 plot of active channels against time
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param scaling scale factor for the data
#' @param sampleHours is the number of hours to plot data for (default is 48)
#' @param sampleIntervalMinutes is the resolution to plot data at
#' @return ggplot2 showing temporal performance
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- SequencingSummaryActiveChannelPlot(seqsum)
#' }
#'
#' @export
SequencingSummaryActiveChannelPlot <- function(seqsum, scaling=1, sampleHours = 48, sampleIntervalMinutes = 60) {

  seqsum$start_time <- seqsum$start_time - min(seqsum$start_time)
  seqsum$start_time <- seqsum$start_time / scaling

  breaks = seq(0, sampleHours*60*60, by=60*sampleIntervalMinutes)
  binass <- findInterval(seqsum$start_time, breaks)

  mergeActiveChannels <- function(interval, binnedAssignments) {
    totalChannels = 0
    if (length(which(binnedAssignments==interval))>0) {
      subset <- seqsum[which(binnedAssignments==interval), ]
      totalChannels = length(unique(subset$channel))
    }
    return(totalChannels)
  }

  binnedTemporalChannels <- data.frame(time=breaks,
                                       channels=unlist(lapply(seq(breaks), mergeActiveChannels, binnedAssignments=binass)
                                       )
  )

  binnedTemporalChannels$time <- binnedTemporalChannels$time / 60 / 60

  activityPlot <- ggplot(binnedTemporalChannels, aes(time)) +
    geom_step(aes(y = channels), size=1, colour = "Steelblue") +
    xlab("Time (hours)") +
    ylab("Number of channels producing data") +
    labs(title="Plot showing number of functional channels against time")
  return(activityPlot)
}








#' present an infographic styled executive summary of sequence_summary.txt content
#'
#' present an infographic styled executive summary of sequence_summary.txt content
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param flowcellId is a label for the plot
#' @return file path to ggplot2 format file
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- SequenceSummaryExecutiveSummary(seqsum)
#' }
#'
#' @export
SequenceSummaryExecutiveSummary <- function(seqsum, flowcellId="undefined") {
  # calculate some basic, but key, metrics

  passedSeqs <- seqsum[which(seqsum$passes_filtering), ]

  readCount <- formatC(nrow(seqsum), big.mark=",")
  totalBases = sum(seqsum$sequence_length_template,na.rm=T)/10^9
  passedBases = sum(passedSeqs$sequence_length_template,na.rm=T)/10^9
  gigabases <- round(totalBases,2)

  # render an info-graphic-like plot for these observations

  infoFile1 <- infoGraphicPlot3(identifier="ExecutiveSummaryValueBoxes",
                                panelA=c(value="flowcell", key=flowcellId, icon="fa-qrcode"),
                                panelB=c(value=readCount, key="Reads produced", icon="fa-filter"),
                                panelC=c(value=gigabases, key="gigabases called", icon="fa-file-text-o"))

  return(infoFile1)

}



#' present an infographic styled basic characteristics plot of sequence_summary.txt content
#'
#' present an infographic styled basic characteristics plot of sequence_summary.txt content
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @return file path to ggplot2 format file
#'
#' @examples
#' \dontrun{
#' seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
#' seqsum <- importSequencingSummary(seqsumFile)
#' plot <- SequenceSummaryBasicInfoPlot(seqsum)
#' }
#'
#' @export
SequenceSummaryBasicInfoPlot <- function(seqsum) {

  passedSeqs <- seqsum[which(seqsum$passes_filtering), ]
  failedSeqs <- seqsum[which(!seqsum$passes_filtering), ]

  passedMeanLength = round(mean(passedSeqs$sequence_length_template), digits = 0)

  N50 <- ncalc(passedSeqs$sequence_length_template, 0.5)

  passedMeanQ = round(mean(passedSeqs$mean_qscore_template), digits = 1)
  failedMeanQ = round(mean(failedSeqs$mean_qscore_template), digits = 1)
  longestRead <- scales::comma_format()(max(passedSeqs$sequence_length_template))

  #N50 length is the length of the shortest contig such that the sum of contigs of equal length or longer is at least 50% of the total length of all contigs

  infoFile2 <- infoGraphicPlot5(identifier="SequenceCharacteristicValueBoxes",
                                panelA=c(value=scales::comma_format()(passedMeanLength), key="Mean Read Length (nt)", icon="fa-bar-chart"),
                                panelB=c(value=scales::comma_format()(N50), key="N50", icon="fa-play"),
                                panelC=c(value=passedMeanQ, key="Mean Read Quality (QV)", icon="fa-area-chart"),
                                panelD=c(value=failedMeanQ, key="Mean Failed QV", icon="fa-bug"),
                                panelE=c(value=longestRead, key="Longest Read", icon="fa-sort"))
  return(infoFile2)
}
