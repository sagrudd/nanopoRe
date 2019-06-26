



#' parse BAM file to calculate summary statistics
#'
#' This method will direct the parsing of a BAM file to extract summary mapping statistics
#'
#' @param bamfilelocation is the location to the BAM file to parse
#' @param window.size is the size of the stepping window to apply (default 100000)
#' @param mc.cores defines the number of cores to use (default = n-1 of available cores)
#' @param force is a logical defining whether the analysis should be recalculated (cached version will be reused if possible)
#' @param level is one of (c("PRIMARY", "SECONDARY", "SUPPLEMENTARY"))
#' @return data.frame of summary observations
#'
#' @examples
#' \dontrun{
#' parseBamFile(file.path("Analysis", "Minimap2", "MyBamFile.bam"))
#' }
#'
#' @export
parseBamFile <- function(bamfilelocation, window.size=100000, mc.cores=min(detectCores()-1, 24), force=FALSE, level="PRIMARY") {
  levels <- c("PRIMARY", "SECONDARY", "SUPPLEMENTARY")
  if (!level %in% levels) {
    warning("LEVEL makes no sense - use PRIMARY|SECONDARY|SUPPLEMENTARY - using PRIMARY INSTEAD")
    level="PRIMARY"
  }
  filenameTag <- ""
  if (level=="SECONDARY") {
    filenameTag <- "secondary_"
  }
  if (level=="SUPPLEMENTARY") {
    filenameTag <- "supplementary_"
  }
  parseBamFileResults <- file.path(getRpath(),
                                   paste0(sub("\\.[^.]*$", "", basename(bamfilelocation)), ".", paste0("aggregated_",filenameTag,as.integer(window.size)), ".Rdata"))
  if (file.exists(parseBamFileResults) & !force) {
    return(readRDS(file=parseBamFileResults))
  }

  parsedBam <- bind_rows(lapply(gtools:::mixedsort(getChromosomeIds()), harvestChromosome, force=force, bamfilelocation=bamfilelocation, window.size=window.size, mc.cores=mc.cores, level=level), .id = "column_label")
  saveRDS(parsedBam, file=parseBamFileResults)
  return(parsedBam)
}



#' parse (indexed) BAM file for summary data for specific chromosome
#'
#' This method will direct the parsing of a specified BAM file for a named chromsome
#'
#' @param chrId is the name of the chromosome of interest
#' @param bamfilelocation is the location to the BAM file to parse
#' @param window.size is the size of the stepping window to apply (default 100000)
#' @param mc.cores defines the number of cores to use (default = n-1 of available cores)
#' @param force is a logical defining whether the analysis should be recalculated (cached version will be reused if possible)
#' @param level is one of (c("PRIMARY", "SECONDARY", "SUPPLEMENTARY"))
#' @return data.frame of summary observations
#'
#' @examples
#' \dontrun{
#' harvestChromosome("1", file.path("Analysis", "Minimap2", "MyBamFile.bam"))
#' }
#'
#' @export
harvestChromosome <- function(chrId, bamfilelocation, window.size=100000, mc.cores=min(detectCores()-1, 24), force=FALSE, level="PRIMARY") {
  levels <- c("PRIMARY", "SECONDARY", "SUPPLEMENTARY")
  if (!level %in% levels) {
    warning("LEVEL makes no sense - use PRIMARY|SECONDARY|SUPPLEMENTARY - using PRIMARY INSTEAD")
    level="PRIMARY"
  }
  filenameTag <- ""
  if (level=="SECONDARY") {
    filenameTag <- "secondary_"
  }
  if (level=="SUPPLEMENTARY") {
    filenameTag <- "supplementary_"
  }

  chromosomeFile <- file.path(getRpath(),
                              paste0(sub("\\.[^.]*$", "", basename(bamfilelocation)), ".", paste0("chrId_",chrId,"_",filenameTag,as.integer(window.size)), ".Rdata"))
  if (file.exists(chromosomeFile) & !force) {
    return(readRDS(file=chromosomeFile))
  }
  dnaStringSetId <- getStringSetId(chrId)

  chromosomeSeq <- getChromosomeSequence(dnaStringSetId)

  coords <- seq(0, length(chromosomeSeq), by=window.size)
  cat(paste("harvest - ", chrId, dnaStringSetId, paste(range(coords), collapse=" "), "\n"))
  # harv <- lapply(coords, harvestBam, dnaStringSetId=dnaStringSetId, chrId=chromosomeId)
  mcharv <- pbmclapply(coords,
                       harvestBam,
                       dnaStringSetId=dnaStringSetId,
                       chrId=chrId,
                       window.size=window.size,
                       bamfilelocation=bamfilelocation,
                       mc.cores=mc.cores,
                       mc.preschedule=FALSE,
                       mc.silent=FALSE,
                       level=level)

  chrData <- as.data.frame(t(as.data.frame(mcharv, stringsAsFactors=FALSE)),
                                             col.names=names(mcharv[[1]]),
                                             row.names=seq_along(length(mcharv)),
                                             stringsAsFactors=FALSE)
  chrData <- fixBamFileColumns(chrData)
  saveRDS(chrData, file=chromosomeFile)
  return(chrData)
}


#' parse (indexed) BAM file for summary data for a given windows on a specific chromosome
#'
#' This method will direct the parsing of a specified BAM file for a named chromsome and qualified interval
#'
#' @param x is a numeric pointer for the interval start position
#' @param dnaStringSetId is the chromosome pointer within the reference genome
#' @param chrId is the name of the chromosome of interest
#' @param window.size is the width of the window in nucleotides
#' @param bamfilelocation is the location to the BAM file to parse
#' @param level is one of (c("PRIMARY", "SECONDARY", "SUPPLEMENTARY"))
#' @return vector of summary observations
#'
#' @examples
#' \dontrun{
#' harvestBam(100000, 1, "1", 100000, file.path("Analysis", "Minimap2", "MyBamFile.bam"))
#' }
#'
#' @export
harvestBam <- function(x, dnaStringSetId, chrId, window.size, bamfilelocation, level="PRIMARY") {
  y = x + window.size
  x <- x + 1

  chromosomeSeq <- getChromosomeSequence(dnaStringSetId)

  if (y > length(chromosomeSeq)) {
    y = length(chromosomeSeq)
  }

  levels <- c("PRIMARY", "SECONDARY", "SUPPLEMENTARY")
  if (!level %in% levels) {
    warning("LEVEL makes no sense - use PRIMARY|SECONDARY|SUPPLEMENTARY - using PRIMARY INSTEAD")
    level="PRIMARY"
  }
  primary <- TRUE
  secondary <- FALSE
  supplementary <- FALSE
  what=c("flag", "strand", "pos", "qwidth", "mapq", "cigar", "qual")
  gccount <- 0
  ncount <- 0
  readq <- 0
  width <- 0

  if (level=="SECONDARY") {
    secondary <- TRUE
    primary <- FALSE
    what=c("strand", "pos", "qwidth", "mapq", "cigar")
  }
  if (level=="SUPPLEMENTARY") {
    supplementary <- TRUE
    primary <- FALSE
    what=c("strand", "pos", "qwidth", "mapq", "cigar")
  }


  params=ScanBamParam(which=GRanges(seqnames = chrId, ranges = IRanges(start = x, end = y)),
                      what=what,
                      flag=scanBamFlag(isSupplementaryAlignment=supplementary, isSecondaryAlignment=secondary),
                      tag=c("NM"))
  SeqCigar <- as.data.frame(scanBam(BamFile(bamfilelocation), param=params)[[1]])

  # to ensure for single counting of reads - select for reads that start within the defined window
  # rather than all reads overlapping the defined window
  onStart <- which(SeqCigar$pos >= x)

  # enumerate the base space which has been mapped from FASTQ
  cigarMapped <- cigarWidthAlongQuerySpace(SeqCigar[onStart, "cigar"], after.soft.clipping=TRUE)

  # count CIGAR derived INS and DEL events
  cigarDEL <- width(cigarRangesAlongPairwiseSpace(SeqCigar[onStart, "cigar"], ops=c("D")))
  cigarINS <- width(cigarRangesAlongPairwiseSpace(SeqCigar[onStart, "cigar"], ops=c("I")))
  cigarDELevents <- unlist(lapply(cigarDEL, length))
  cigarINSevents <- unlist(lapply(cigarINS, length))
  cigarDELbases <- unlist(lapply(cigarDEL, sum))
  cigarINSbases <- unlist(lapply(cigarINS, sum))

  if (primary) {
    letterFreq <- letterFrequency(subseq(chromosomeSeq, x, y), c("A", "C", "G", "T", "N"))
    width <- sum(letterFreq)

    # parse out the mean per-read base calling qvalues
    #alphabetScore <- alphabetScore(FastqQuality(SeqCigar[onStart, "qual"])) / SeqCigar[onStart, "qwidth"]
    alphabetScore <- unlist(lapply(as.character(SeqCigar[onStart, "qual"]), qualToMeanQ))

    gccount <- as.integer(letterFreq['G'] + letterFreq['C'])
    ncount <- as.integer(letterFreq['N'])
    readq <- round(phredmean(alphabetScore), digits=2)
  }

  # process the depth of coverage stuff ...
  tb <- c(0,0,0,0,0)
  meanCov <- 0
  if (nrow(SeqCigar)>0) {
    cov <- GAlignments(seqnames=rep("map", nrow(SeqCigar)), pos=SeqCigar$pos, strand=SeqCigar$strand, cigar=as.character(SeqCigar$cigar))
    covTable <- as.data.frame(coverage(cov, shift=-x, width=window.size))
    tb <- quantile(covTable$value)
    meanCov <- mean(covTable$value)
  }
  names(tb) <- c("iqr.0", "iqr.25", "iqr.50", "iqr.75", "iqr.100")

  res <- c(chrId=chrId,
           startpos=x,
           readStarts=length(onStart),                        # validated as correct with Chr20 dataset (samtools as ground truth)
           plusStrand=length(which(SeqCigar[onStart, "strand"]=="+")),
           basesReadsStarted=sum(SeqCigar[onStart,"qwidth"]), # validated as correct with Chr20 dataset
           width=width,                                       # the width of the window considered
           gccount=gccount,
           ncount=ncount,
           mismatches=sum(SeqCigar[onStart,"NM"]),            # the sum of total mapping edit distance
           cigarMapped=sum(cigarMapped),                      # the sum of read bases that are mapped to reference genome, starting in this window, after softclipping
           cigarInsertionBases=sum(cigarINSbases),            # the sum of inserted bases against the reference genome (pairwise presentation)
           cigarDeletionBases=sum(cigarDELbases),             # the sum of inserted bases against the reference genome (pairwise presentation)
           cigarInsertionEvents=sum(cigarINSevents),
           cigarDeletionEvents=sum(cigarDELevents),
           mapq=round(phredmean(SeqCigar[onStart,"mapq"]), digits=2),# mean mapping quality for reads starting in this window
           readQ=readq,        # this is the per-read mean Q value averaged across onStart
           readLen=round(mean(mean(SeqCigar[onStart, "qwidth"])),digits=2),  # average sequence length for onStart reads
           unlist(tb),                                                 # the depth of coverage summary
           meanCov=round(meanCov, digits=2),
           flagCount=length(unique(SeqCigar$flag))
  )
  return(res)
}


#' extract unmapped read quality information from provided BAM file
#'
#' This method will parse the unmapped reads from a BAM file and return a data frame containing
#' readID and quality information
#'
#' @param bamfilelocation is the location to the BAM file to parse
#' @param force is a boolean defining whether an existing result file should be overwritten
#' @return vector of summary observations
#'
#' @examples
#' \dontrun{
#' harvestUnmappedReads(file.path("Analysis", "Minimap2", "MyBamFile.bam"))
#' }
#'
#' @export
harvestUnmappedReads <- function(bamfilelocation, force=FALSE) {
  # this is a slow method since the whole file needs to be parsed ... - the parsing logic is not
  # perfect and there is considerable system IO ...

  # alternative samtools view -@ 4 -f 0x4 -O SAM ./Analysis/GM24385.15/alignment/GM24385.15_minimap2.bam | awk '{{print $11}}'

  unmappedReads <- file.path(getRpath(),
                             paste0(sub("\\.[^.]*$", "", basename(bamfilelocation)), "_UnmappedReads", ".Rdata"))
  if (file.exists(unmappedReads) & !force) {
    return(readRDS(file=unmappedReads))
  }

  params=ScanBamParam(what=c("qname", "flag", "qual"),
                      flag=scanBamFlag(isUnmappedQuery=TRUE))
  unmappedReadData <- as.data.frame(scanBam(BamFile(bamfilelocation), param=params)[[1]])
  saveRDS(unmappedReadData, file=unmappedReads)
  return(unmappedReadData)

}






fixBamFileColumns <- function(parsedBamFile) {
  if (!is.integer(parsedBamFile$startpos)) { parsedBamFile$startpos <- as.integer(parsedBamFile$startpos) }
  if (!is.integer(parsedBamFile$readStarts)) { parsedBamFile$readStarts <- as.integer(parsedBamFile$readStarts) }
  if (!is.integer(parsedBamFile$plusStrand)) { parsedBamFile$plusStrand <- as.integer(parsedBamFile$plusStrand) }
  if (!is.integer(parsedBamFile$basesReadsStarted)) { parsedBamFile$basesReadsStarted <- as.integer(parsedBamFile$basesReadsStarted) }
  if (!is.integer(parsedBamFile$width)) { parsedBamFile$width <- as.integer(parsedBamFile$width) }
  if (!is.integer(parsedBamFile$gccount)) { parsedBamFile$gccount <- as.integer(parsedBamFile$gccount) }
  if (!is.integer(parsedBamFile$ncount)) { parsedBamFile$ncount <- as.integer(parsedBamFile$ncount) }
  if (!is.integer(parsedBamFile$mismatches)) { parsedBamFile$mismatches  <- as.integer(parsedBamFile$mismatches) }
  if (!is.integer(parsedBamFile$cigarMapped)) { parsedBamFile$cigarMapped  <- as.integer(parsedBamFile$cigarMapped) }
  if (!is.integer(parsedBamFile$cigarInsertionBases)) { parsedBamFile$cigarInsertionBases  <- as.integer(parsedBamFile$cigarInsertionBases) }
  if (!is.integer(parsedBamFile$cigarDeletionBases)) { parsedBamFile$cigarDeletionBases  <- as.integer(parsedBamFile$cigarDeletionBases) }
  if (!is.integer(parsedBamFile$cigarInsertionEvents)) { parsedBamFile$cigarInsertionEvents  <- as.integer(parsedBamFile$cigarInsertionEvents) }
  if (!is.integer(parsedBamFile$cigarDeletionEvents)) { parsedBamFile$cigarDeletionEvents  <- as.integer(parsedBamFile$cigarDeletionEvents) }
  if (!is.integer(parsedBamFile$iqr.0)) { parsedBamFile$iqr.0  <- as.integer(parsedBamFile$iqr.0) }
  if (!is.integer(parsedBamFile$iqr.25)) { parsedBamFile$iqr.25  <- as.integer(parsedBamFile$iqr.25) }
  if (!is.integer(parsedBamFile$iqr.50)) { parsedBamFile$iqr.50  <- as.integer(parsedBamFile$iqr.50) }
  if (!is.integer(parsedBamFile$iqr.75)) { parsedBamFile$iqr.75  <- as.integer(parsedBamFile$iqr.75) }
  if (!is.integer(parsedBamFile$iqr.100)) { parsedBamFile$iqr.100  <- as.integer(parsedBamFile$iqr.100) }

  if (!is.numeric(parsedBamFile$mapq)) { parsedBamFile$mapq  <- as.numeric(parsedBamFile$mapq) }
  if (!is.numeric(parsedBamFile$readQ)) { parsedBamFile$readQ  <- as.numeric(parsedBamFile$readQ) }
  if (!is.numeric(parsedBamFile$readLen)) { parsedBamFile$readLen  <- as.numeric(parsedBamFile$readLen) }
  if (!is.numeric(parsedBamFile$meanCov)) { parsedBamFile$meanCov  <- as.numeric(parsedBamFile$meanCov) }

  return(parsedBamFile)
}



