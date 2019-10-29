

#' @importFrom rhdf5 h5readAttributes
#' @importFrom xts .parseISO8601
#' @importFrom dplyr last
extractH5Modstemplate <- function(read_id, h5content, fast5file) {
    sample <- h5content[grep(read_id, h5content$group),]
    # get the BaseModProbs ...
    modBaseProbs <- which(sample$name=="ModBaseProbs")
    if (length(modBaseProbs)==0) {
        stop(paste0("No Base Modification table in Fast5 file element [",sample,"]"))
    } else if (length(modBaseProbs)==1) {
        return(basename(dirname(sample[modBaseProbs,"group"])))
    } else {
        # there are multiple methylation calls in this dataset - choose the most recent
        # https://en.wikipedia.org/wiki/ISO_8601#Combined_date_and_time_representations
        times <- lapply(modBaseProbs, function(x) {
            xts::.parseISO8601(
                gsub("Z$", "", rhdf5::h5readAttributes(fast5file, dirname(sample[x,"group"]))$time_stamp), tz="UTC")$first.time})
        pos <- modBaseProbs[dplyr::last(order(unlist(times)))]
        return(basename(dirname(sample[pos,"group"])))
    }
}




#' prepares a data.frame of per-read methylC bases from provided fast5 file
#'
#' This method will parse a fast5 file to extract the 5mC bases classified as 5mC with a
#' probability >= the provided threshold_5mc variable
#'
#' @param fast5file is the fast5 file to parse
#' @param threshold_5mc is the cutoff to apply to values returned
#' @param fast5files is a vector of files being processed (for pretty logging)
#' @param mc.cores the number of cores to use on a multi-processor Linux system
#' @param force boolean as to whether to force the calculation
#' @param ... beyond
#'
#' @return data.frame of coordinates and sequence names from fast5
#'
#' @export
extract5mC <- function(fast5file, threshold_5mc = 0.85, fast5files=NULL, force = FALSE, ...) {

    mod_data <- extractModifiedBasesFromFast5(fast5file, fast5files, mc.cores)
    mod_data <- mod_data[which(mod_data$prob_5mC>=threshold_5mc),]

    return(invisible(mod_data))
}





#' prepares a data.frame of read mapping context for methylation parsing
#'
#' This method will parses the specified BAM file to extract CIGAR and other coordinates for
#' the purpose of selecting sequence bases that are base modified
#'
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBam
#' @importFrom Rsamtools bamFlagAsBitMatrix
#' @param bamFile is the bam file to parse
#' @param chrId is the identifier for the chromosome to parse
#' @param force is logical specifying whether the load should be forced (recalculated)
#'
#' @return data.frame
#'
#' @export
loadMethylationBamFile <- function(bamFile, chrId, force=FALSE) {
    print(paste0("parsing bam [",bamFile,"@",chrId,"]"))

    bamParseResults <- file.path(getRpath(), paste0(digest::digest(bamFile, algo="md5", file = FALSE), ".chr.", chrId, ".Rdata"))
    if (file.exists(bamParseResults) && !force) {
        bamChunk <- readRDS(file = bamParseResults)
        return(invisible(bamChunk))
    }

    what = c("qname", "flag", "rname", "strand", "pos", "qwidth", "cigar")
    param = Rsamtools::ScanBamParam(which = GenomicRanges::GRanges(seqnames = chrId, ranges = IRanges::IRanges(start = 1, end = as.numeric(getSeqLengths(chrId)))),
                                    what = what)

    bam <- open(Rsamtools::BamFile(bamFile))
    reads = Rsamtools::scanBam(bam, param = param)[[1L]]

    close(bam)

    # dispose of the secondary, supplementary and unmapped reads
    bamChunk <- as.data.frame(reads)
    flagMatrix <- as.data.frame(Rsamtools::bamFlagAsBitMatrix(as.integer(bamChunk$flag)))

    secondary <- which(flagMatrix$isSecondaryAlignment==1)
    supplemen <- which(flagMatrix$isSupplementaryAlignment==1)
    unmapped <- which(flagMatrix$isUnmappedQuery==1)

    drop <- unique(append(append(secondary, supplemen), unmapped))
    bamChunk <- bamChunk[-drop, ]

    saveRDS(bamChunk, file = bamParseResults)

    return(invisible(bamChunk))
}




#' prepares a data.frame of read mapping context for methylation parsing
#'
#' This method will parses the specified BAM file to extract CIGAR and other coordinates for
#' the purpose of selecting sequence bases that are base modified
#'
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBam
#' @importFrom Rsamtools bamFlagAsBitMatrix
#' @param x is a numeric pointer
#' @param filteredChunk corresponds to parsed BAM file data.frame
#' @param modifications_df corresponds to parsed methylation data.frame
#'
#' @return data.frame
#'
#' @export
extractMethylatedBases <- function(x, filteredChunk, modifications_df) {
    mapData <- filteredChunk[x,]
    variants <- which(modifications_df$read_id == mapData$qname)
    varData <- modifications_df[variants,]

    if (mapData$strand == '+') {
        referenceCoordinates <- cigarQ2R(varData$position, mapData$cigar)  + mapData$pos
        op <- cigarQ2Op(varData$position, mapData$cigar)
    } else {
        referenceCoordinates <- cigarQ2R((mapData$qwidth - varData$position - 1), mapData$cigar)  + mapData$pos
        op <- cigarQ2Op((mapData$qwidth - varData$position - 1), mapData$cigar)
    }

    # there is a sneaky exception case where NA values creep in ...
    if (NA %in% referenceCoordinates) {
        op <- op[-which(is.na(referenceCoordinates))]
        referenceCoordinates <- referenceCoordinates[-which(is.na(referenceCoordinates))]
    }

    pdata <- NA

    if (length(referenceCoordinates)>0) {
        bases <- getChromosomeSequence(getStringSetId(mapData$rname))[referenceCoordinates]
        if (mapData$strand == '-') {
            bases <- Biostrings::complement(bases)
        }
        rnt <- strsplit(as.character(bases), "")[[1]]
        pdata <- data.frame(chr=as.character(mapData$rname),
                            pos=referenceCoordinates,
                            fwd=as.integer(mapData$strand == '+'),
                            rev=as.integer(mapData$strand == '-'),
                            op=as.character(op),
                            A=as.integer(rnt=='A'),
                            C=as.integer(rnt=='C'),
                            G=as.integer(rnt=='G'),
                            T=as.integer(rnt=='T'),
                            stringsAsFactors = FALSE)
        pdata <- pdata[which(pdata$op == "M"), ]
    }

    return(invisible(pdata))
}




extractModifiedBasesFromFast5 <- function(fast5file, fast5files=NULL, mc.cores=(parallel::detectCores()-1)) {
    if (is.null(fast5files)) {
        cat(paste0("extracting 5mC probabilities from [",fast5file,"]\n"))
    } else {
        cat(paste0("[",which(fast5files==fast5file)," of ",length(fast5files),"] 5mC extract [",fast5file,"]\n"))
    }

    baseModResults <- file.path(getRpath(), paste0(digest::digest(fast5file, algo="md5", file = FALSE), ".basemods", ".Rdata"))
    if (file.exists(baseModResults) && !force) {
        return(invisible(readRDS(file = baseModResults)))
    }
    h5content <- rhdf5::h5ls(fast5file)
    read_ids <- h5content[which(h5content$group == "/"),"name"]

    # use the first read_id to parse out dates of
    read_id <- read_ids[1]
    methtemplate <- extractH5Modstemplate(read_id, h5content=h5content, fast5file=fast5file)

    extract5mCByRead <- function(read_id, fast5file, template="Basecall_1D_000") {
        targetFastq <- paste0("/",read_id,"/Analyses/",template,"/BaseCalled_template/Fastq")
        targetMods <- paste0("/",read_id,"/Analyses/",template,"/BaseCalled_template/ModBaseProbs")

        read_id <- gsub("(^/read_)|(/Analyses.+)", "", targetFastq)

        fastq <- unlist(strsplit(unlist(strsplit(rhdf5::h5read(fast5file, targetFastq), "\n"))[2], ""))

        basemodprob <- data.frame(t(rhdf5::h5read(fast5file, targetMods)))[,4]/255

        mods_df <- data.frame(read_id=read_id, position=seq(length(fastq)), nucleotide=fastq, prob_5mC=basemodprob, stringsAsFactors=FALSE)
        mods_df <- mods_df[which(mods_df$nucleotide=="C"),]
        return(invisible(mods_df))
    }

    mod_data <- dplyr::bind_rows(pbmcapply::pbmclapply(read_ids, extract5mCByRead, fast5file=fast5file, template=methtemplate, mc.cores=mc.cores))
    saveRDS(mod_data, file = baseModResults)
    return(mod_data)
}



#' prepares a data.frame of per-read methylC probability counts
#'
#' This method will parse a fast5 file to the distribution of probabilities for C residues
#' actually being 5mC residues
#'
#' @param fast5file is the fast5 file to parse
#' @param fast5files is a vector of files being processed (for pretty logging)
#' @param force boolean as to whether to force the calculation
#' @param ... beyond
#'
#' @return data.frame of coordinates and sequence names from fast5
#'
#' @export
extract5mCProbabilities <- function(fast5file, fast5files=NULL, force = FALSE, ...) {
    mod_data <- extractModifiedBasesFromFast5(fast5file, fast5files)
    tab <- as.integer(table(factor(round(mod_data$prob_5mC * 100), levels = seq(0, 100))))
    return(invisible(tab))
}
