


#' extract content from a BAM file
#'
#' This method will extract a block of BAM content from a given BAM file
#'
#' @param bamFile is the location to the BAM file to parse
#' @param yieldSize is the number of mapping observations to extract
#' @return list of summary observations
#'
#' @examples
#' # define the path to a BAM file
#' demoBam <- system.file("extdata",
#'     "Ecoli_zymo_R10_filt_subs.bam",
#'     package = "nanopoRe")
#' bamdata <- testBam(demoBam, 10)
#'
#' @export
testBam <- function(bamFile, yieldSize = 100L) {
    bam = open(BamFile(bamFile, yieldSize = yieldSize))
    what = c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "qual", "cigar")
    param = ScanBamParam(what = what, tag = c("NM", "MD"))
    reads = scanBam(bam, param = param)[[1L]]
    close(bam)
    return(reads)
}



#' summarise mapping observations from a BAM file by chromosome
#'
#' This method will parse a BAM file and summarise the mapping observations in a
#' way that can be used for preparation of figures and confidence plots by
#' the nanopoRe package
#'
#' @usage bamSummariseByChr(chrId,
#'     bamFile,
#'     force = FALSE,
#'     blockSize = 50000L,
#'     index=NULL
#' )
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBam
#' @param chrId is the chromosome to parse
#' @param bamFile is the location to the BAM file to parse
#' @param force logical value describing whether the analysis should be force recalculated
#' @param blockSize the number of reads to process at the time as iterating through
#' @param index path to the BAI index file - should be automatic in most cases ...
#' @return data.frame of per read summary observations
#'
#' @examples
#' # load reference genome
#' init()
#' referenceFasta <- system.file("extdata",
#'     "Escherichia_coli_complete_genome.fasta",
#'     package = "nanopoRe")
#' setReferenceGenome(referenceFasta)
#' # load reference BAM
#' demoBam <- system.file("extdata",
#'     "Ecoli_zymo_R10_filt_subs.bam",
#'     package = "nanopoRe")
#' # Rsamtools::BamFile has problems in picking up packaged .bam.bai
#' # we will specify it here rather then automatic detection
#' demoBamIdx <- system.file("extdata",
#'     "Ecoli_zymo_R10_filt_subs.bam",
#'     package = "nanopoRe")
#' bamSummariseByChr(chrId=getChromosomeIds()[1],
#'     bamFile=demoBam,
#'     blockSize=100L,
#'     index=demoBamIdx)
#'
#' @export
bamSummariseByChr <- function(chrId, bamFile, force = FALSE, blockSize = 50000L, index=NULL) {
    bamSummaryResults <- file.path(getRpath(), paste0(sub("\\.[^.]*$", "", basename(bamFile)), ".bamMetrics.chr",
        chrId, ".Rdata"))
    message(paste0("targetRdata ", bamSummaryResults, "\n"))
    if (file.exists(bamSummaryResults) & !force) {
        return(readRDS(file = bamSummaryResults))
    }
    count <- 0
    bamInfo <- data.frame()
    bam <- NULL
    if (!is.null(index)) {
        bam = open(BamFile(bamFile, yieldSize = blockSize, index=index))
    } else {
        bam = open(BamFile(bamFile, yieldSize = blockSize))
    }
    what = c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "qual", "cigar")
    param = ScanBamParam(which = GRanges(seqnames = chrId, ranges = IRanges(start = 1, end = as.numeric(getSeqLengths(chrId)))),
        what = what, tag = c("NM", "MD"))
    repeat {
        reads = scanBam(bam, param = param)[[1L]]
        if (length(reads$qname) == 0L)
            break
        count = count + length(reads$qname)
        cat(paste0(count, "\n"))
        bamInfo <- rbind(bamInfo, processBamChunk(reads))
    }
    close(bam)
    saveRDS(bamInfo, file = bamSummaryResults)
    return(bamInfo)
}



#' summarise mapping observations from a BAM file using parallel chr-by-chr
#'
#' This method will parse a BAM file and summarise the mapping observations in a
#' way that can be used for preparation of figures and confidence plots by
#' the nanopoRe package
#'
#' @importFrom parallel detectCores
#' @importFrom pbmcapply pbmclapply
#' @usage parallelBamSummarise(bamFile, force = FALSE,
#'     mc.cores = min(parallel::detectCores() - 1, 24))
#' @param bamFile is the location to the BAM file to parse
#' @param force logical value describing whether the analysis should be force recalculated
#' @param mc.cores number of threads to use for the process
#' @return data.frame of per read summary observations
#'
#' @examples
#' \dontrun{
#' parallelBamSummarise(file.path('Analysis', 'Minimap2', 'MyBamFile.bam'))
#' }
#'
#' @export
parallelBamSummarise <- function(bamFile, force = FALSE, mc.cores = min(parallel::detectCores() - 1,
    24)) {
    bamSummaryResults <- file.path(getRpath(), paste0(sub("\\.[^.]*$", "", basename(bamFile)), ".bamMetrics",
        ".Rdata"))
    message(paste0("targetRdata ", bamSummaryResults, "\n"))
    if (file.exists(bamSummaryResults) & !force) {
        return(readRDS(file = bamSummaryResults))
    }
    mcharv <- pbmclapply(getChromosomeIds(), bamSummariseByChr, bamFile = bamFile, force = force, mc.cores = mc.cores,
        mc.preschedule = FALSE, mc.silent = FALSE)
    result <- data.frame()
    for (chr in getChromosomeIds()) {
        result <- rbind(result, bamSummariseByChr(chr, bamFile))
    }
    saveRDS(result, file = bamSummaryResults)
    return(result)
}



#' summarise mapping observations from a BAM file
#'
#' This method will parse a BAM file and summarise the mapping observations in a
#' way that can be used for preparation of figures and confidence plots by
#' the nanopoRe package
#'
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBam
#' @param bamFile is the location to the BAM file to parse
#' @param force logical value describing whether the analysis should be force recalculated
#' @param blockSize the number of reads to process at the time as iterating through
#' @return data.frame of per read summary observations
#'
#' @examples
#' demoBam <- system.file("extdata",
#'     "Ecoli_zymo_R10_filt_subs.bam",
#'     package = "nanopoRe")
#' bamSummary <- bamSummarise(demoBam, force=FALSE, blockSize=1000L)
#'
#'
#' @export
bamSummarise <- function(bamFile, force = FALSE, blockSize = 50000L) {
    bamSummaryResults <- file.path(getRpath(), paste0(sub("\\.[^.]*$", "", basename(bamFile)), ".bamMetrics",
        ".Rdata"))
    message(paste0("targetRdata ", bamSummaryResults, "\n"))
    if (file.exists(bamSummaryResults) & !force) {
        return(getCachedFileObject("BamTargetData", bamSummaryResults))
    }
    count <- 0
    bamInfo <- data.frame()
    bam = open(BamFile(bamFile, yieldSize = blockSize))
    what = c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "qual", "cigar")
    param = ScanBamParam(what = what, tag = c("NM", "MD"))
    repeat {
        reads = scanBam(bam, param = param)[[1L]]
        if (length(reads$qname) == 0L)
            break
        count = count + length(reads$qname)
        cat(paste0(count, "\n"))
        bamInfo <- rbind(bamInfo, processBamChunk(reads))
    }
    close(bam)
    saveRDS(bamInfo, file = bamSummaryResults)
    return(bamInfo)
}


#' @import Rsamtools
#' @importFrom GenomicAlignments cigarRangesAlongPairwiseSpace
#' @importFrom GenomicAlignments cigarRangesAlongQuerySpace
#' @importFrom GenomicAlignments cigarRangesAlongReferenceSpace
#' @importFrom IRanges width
processBamChunk <- function(bamChunk) {
    bamChunk <- as.data.frame(bamChunk)
    # separate out unmapped reads
    flagMatrix <- as.data.frame(bamFlagAsBitMatrix(as.integer(bamChunk$flag)))
    unmappedChunk <- bamChunk[which(flagMatrix$isUnmappedQuery == 1), ]
    bamChunk <- bamChunk[which(flagMatrix$isUnmappedQuery == 0), ]

    # simplify the read flags ...
    readFlag <- factor(rep("Primary", length(bamChunk$flag)), levels = c("Primary", "Secondary", "Supplementary",
        "Unmapped"))
    flagMatrix <- as.data.frame(bamFlagAsBitMatrix(as.integer(bamChunk$flag)))
    readFlag[which(flagMatrix$isSecondaryAlignment == 1)] <- "Secondary"
    readFlag[which(flagMatrix$isSupplementaryAlignment == 1)] <- "Supplementary"
    readFlag[which(flagMatrix$isUnmappedQuery == 1)] <- "Unmapped"
    # calculate aln_len for alignment
    ref_aln_len <- unlist(width(cigarRangesAlongReferenceSpace(bamChunk$cigar, reduce.ranges = TRUE)))
    endpos <- bamChunk$pos + ref_aln_len
    # mean readq
    readq <- unlist(lapply(as.character(bamChunk$qual), qualToMeanQ))

    # get the actual query bases mapped ...
    q_aln_len <- unlist(IRanges::width(cigarRangesAlongQuerySpace(bamChunk$cigar, after.soft.clipping = TRUE,
        reduce.ranges = TRUE)))
    # q_ins_len <- sum(IRanges::width(cigarRangesAlongQuerySpace(bamChunk$cigar, drop.empty.ranges =
    # TRUE, after.soft.clipping=TRUE, ops=c('I'), with.ops = TRUE)))
    q_ins_len <- sum(IRanges::width(cigarRangesAlongPairwiseSpace(bamChunk$cigar, ops = c("I"))))
    q_del_len <- sum(IRanges::width(cigarRangesAlongPairwiseSpace(bamChunk$cigar, ops = c("D"))))
    q_match_len <- sum(IRanges::width(cigarRangesAlongPairwiseSpace(bamChunk$cigar, ops = c("M"))))
    q_mm_len <- bamChunk$tag.NM - q_ins_len - q_del_len
    coverage = q_aln_len/bamChunk$qwidth
    accuracy = (q_match_len - q_mm_len)/(q_match_len + q_ins_len + q_del_len)
    identity = (q_match_len - q_mm_len)/(q_match_len)

    parsed <- data.frame(qname = bamChunk$qname, readFlag = readFlag, rname = bamChunk$rname, strand = bamChunk$strand,
        start = bamChunk$pos, qwidth = bamChunk$qwidth, end = endpos, mapq = bamChunk$mapq, readq = readq,
        coverage = coverage, accuracy = accuracy, identity = identity)

    if (nrow(unmappedChunk) > 0) {
        unmapParse <- data.frame(qname = unmappedChunk$qname, readFlag = factor(rep("Unmapped", length(unmappedChunk$flag)),
            levels = c("Primary", "Secondary", "Supplementary", "Unmapped")), rname = NA, strand = factor(rep("*",
            length(unmappedChunk$strand)), levels = c("+", "-", "*")), start = NA, end = NA, qwidth = nchar(unmappedChunk$qual),
            mapq = NA, readq = unlist(lapply(as.character(unmappedChunk$qual), qualToMeanQ)), coverage = NA,
            accuracy = NA, identity = NA)
        parsed <- rbind(parsed, unmapParse)
    }
    return(parsed)
}


#' get depth of coverage information from across the genome
#'
#' This method will return a tiled Granges object containing mean depth of coverage information
#'
#' @usage bamSummaryToCoverage(bamFile,
#'     tilewidth = 1e+05,
#'     blocksize = 10000,
#'     flag = 'Primary'
#' )
#' @param bamFile - path to the bamFile to use
#' @param tilewidth - the size of the window to use for the tiling
#' @param blocksize to use for parsing the BAM file
#' @param flag the mapping type to filter reads for (Primary/Secondary/Supplementary)
#' @return GRanges object with mean depth of coverage data in binned_cov field
#'
#' @examples
#' demoBam <- system.file("extdata", "Ecoli_zymo_R10_filt_subs.bam", package = "nanopoRe")
#' bamGR <- bamSummaryToCoverage(demoBam)
#'
#' @export
bamSummaryToCoverage <- function(bamFile, tilewidth = 1e+05, blocksize = 10000, flag = "Primary") {
    bamSummary <- bamSummarise(bamFile, blockSize = 10000)
    primary <- bamSummary[which(bamSummary$readFlag == flag), ]

    # depending on the genome used there may be a load of warnings here this is likely due to reads
    # mapping beyond segment boundaries - warnings are masked here since they are expected
    suppressWarnings(grdata <- GRanges(seqnames = primary$rname, ranges = IRanges(start = primary$start,
        end = primary$end), strand = primary$strand, seqlengths = getSeqLengths(levels(primary$rname))))
    mapCoverage <- coverage(grdata)
    bins <- tileGenome(seqlengths(grdata), tilewidth=tilewidth, cut.last.tile.in.chrom = TRUE)
    return(binnedAverage(bins, mapCoverage, "binned_cov"))
}


