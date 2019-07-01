


#' extract content from a BAM file
#'
#' This method will extract a block of BAM content from a given BAM file
#'
#' @param bamFile is the location to the BAM file to parse
#' @return list of summary observations
#'
#' @examples
#' \dontrun{
#' testBam(file.path("Analysis", "Minimap2", "MyBamFile.bam"))
#' }
#'
#' @export
testBam <- function(bamFile) {
  bam = open(BamFile(bamFile, yieldSize = 10000L))
  what=c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "qual", "cigar")
  param=ScanBamParam(what=what, tag=c("NM", "MD"))
  reads = scanBam(bam, param = param)[[1L]]
  close(bam)
  return(reads)
}


#' summarise mapping observations from a BAM file
#'
#' This method will parse a BAM file and summarise the mapping observations in a
#' way that can be used for preparation of figures and confidence plots by
#' the nanopoRe package
#'
#' @param bamFile is the location to the BAM file to parse
#' @param force logical value describing whether the analysis should be force recalculated
#' @param blockSize the number of reads to process at the time as iterating through
#' @return data.frame of per read summary observations
#'
#' @examples
#' \dontrun{
#' bamSummarise(file.path("Analysis", "Minimap2", "MyBamFile.bam"), force=FALSE, blockSize=1000L)
#' }
#'
#' @export
bamSummarise <- function(bamFile, force=FALSE, blockSize=50000L) {
  bamSummaryResults <- file.path(getRpath(),
                                 paste0(sub("\\.[^.]*$", "", basename(bamFile)), ".bamMetrics", ".Rdata"))
  message(paste0("targetRdata ",bamSummaryResults,"\n"))
  if (file.exists(bamSummaryResults) & !force) {
    return(readRDS(file=bamSummaryResults))
  }
  count <- 0
  bamInfo <- data.frame()
  bam = open(BamFile(bamFile, yieldSize = blockSize))
  what=c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "qual", "cigar")
  param=ScanBamParam(what=what, tag=c("NM", "MD"))
  repeat {
    reads = scanBam(bam, param = param)[[1L]]
    if (length(reads$qname) == 0L) break
    count = count + length(reads$qname)
    cat(paste0(count, "\n"))
    bamInfo <- rbind(bamInfo, processBamChunk(reads))
  }
  close(bam)
  cat(paste0(count, "\n"))
  saveRDS(bamInfo, file=bamSummaryResults)
  return(bamInfo)
}



processBamChunk <- function(bamChunk) {
  bamChunk <- as.data.frame(bamChunk)

  # separate out unmapped reads
  flagMatrix <- as.data.frame(bamFlagAsBitMatrix(as.integer(bamChunk$flag)))
  unmappedChunk <- bamChunk[which(flagMatrix$isUnmappedQuery==1),]
  bamChunk <-  bamChunk[which(flagMatrix$isUnmappedQuery==0),]

  # simplify the read flags ...
  readFlag <- factor(rep("Primary",length(bamChunk$flag)), levels=c("Primary", "Secondary", "Supplementary", "Unmapped"))
  flagMatrix <- as.data.frame(bamFlagAsBitMatrix(as.integer(bamChunk$flag)))
  readFlag[which(flagMatrix$isSecondaryAlignment==1)] <- "Secondary"
  readFlag[which(flagMatrix$isSupplementaryAlignment==1)] <- "Supplementary"
  readFlag[which(flagMatrix$isUnmappedQuery==1)] <- "Unmapped"

  # calculate aln_len for alignment
  ref_aln_len <- unlist(width(cigarRangesAlongReferenceSpace(bamChunk$cigar, reduce.ranges = TRUE)))
  endpos <- bamChunk$pos + ref_aln_len

  # mean readq
  readq <- unlist(lapply(as.character(bamChunk$qual), qualToMeanQ))

  # get the actual query bases mapped ...
  q_aln_len <- unlist(width(cigarRangesAlongQuerySpace(bamChunk$cigar, after.soft.clipping=TRUE, reduce.ranges = TRUE)))
  #q_ins_len <- sum(width(cigarRangesAlongQuerySpace(bamChunk$cigar, drop.empty.ranges = TRUE, after.soft.clipping=TRUE, ops=c("I"), with.ops = TRUE)))
  q_ins_len <- sum(width(cigarRangesAlongPairwiseSpace(bamChunk$cigar, ops=c("I"))))
  q_del_len <- sum(width(cigarRangesAlongPairwiseSpace(bamChunk$cigar, ops=c("D"))))
  q_match_len <- sum(width(cigarRangesAlongPairwiseSpace(bamChunk$cigar, ops=c("M"))))
  q_mm_len <- bamChunk$tag.NM - q_ins_len - q_del_len
  coverage = q_aln_len / bamChunk$qwidth
  accuracy = (q_match_len - q_mm_len) / (q_match_len)
  identity = (q_match_len - q_mm_len) / (q_match_len + q_ins_len)

  parsed <- data.frame(
    qname=bamChunk$qname,
    readFlag=readFlag,
    rname=bamChunk$rname,
    strand=bamChunk$strand,
    start=bamChunk$pos,
    qwidth=bamChunk$qwidth,
    end=endpos,
    mapq=bamChunk$mapq,
    readq=readq,
    coverage=coverage,
    accuracy=accuracy,
    identity=identity
  )

  if (nrow(unmappedChunk) > 0) {
    unmapParse <- data.frame(
      qname=unmappedChunk$qname,
      readFlag=factor(rep("Unmapped",length(unmappedChunk$flag)), levels=c("Primary", "Secondary", "Supplementary", "Unmapped")),
      rname=NA,
      strand=factor(rep("*",length(unmappedChunk$strand)), levels=c("+", "-", "*")),
      start=NA,
      end=NA,
      qwidth=nchar(unmappedChunk$qual),
      mapq=NA,
      readq=unlist(lapply(as.character(unmappedChunk$qual), qualToMeanQ)),
      coverage=NA,
      accuracy=NA,
      identity=NA
    )
    parsed <- rbind(parsed, unmapParse)

  }
  return(parsed)
}


#' get depth of coverage information from across the genome
#'
#' This method will return a tiled Granges object containing mean depth of coverage information
#'
#' @param bamFile - path to the bamFile to use
#' @param tilewidth - the size of the window to use for the tiling
#' @param blocksize to use for parsing the BAM file
#' @param flag the mapping type to filter reads for (Primary/Secondary/Supplementary)
#' @return GRanges object with mean depth of coverage data in binned_cov field
#'
#' @examples
#' \dontrun{
#' bamSummaryToCoverage(file.path("Analysis", "Minimap2", "MyBamFile.bam"))
#' }
#'
#' @export
bamSummaryToCoverage <- function(bamFile, tilewidth=100000, blocksize=10000, flag="Primary") {
  bamSummary <- bamSummarise(bamFile, blockSize=10000)
  primary <- bamSummary[which(bamSummary$readFlag==flag),]

  grdata <- GRanges(seqnames=primary$rname,
                    ranges=IRanges(start=primary$start, end=primary$end),
                    strand=primary$strand,
                    seqlengths=getSeqLengths(levels(primary$rname)))
  mapCoverage <- coverage(grdata)
  bins <- tileGenome(seqlengths(grdata), tilewidth=100000, cut.last.tile.in.chrom=TRUE)
  return(binnedAverage(bins, mapCoverage, "binned_cov"))
}


