

#' perform an R-based enrichment analysis for a target region of the genome
#'
#' This method is associated with the cas9-mediated enrichment of genomic
#' target regions ...
#'
#' @param force defines whether a reanalysis of the data should be forced
#' if results already existing
#' @param ... - who knows - I'm using this for mc.cores
#' @return TRUE or FALSE depending
#'
#' @examples
#' yamlFile <- system.file("extdata", "cas9_demo.yaml", package = "nanopoRe")
#' sourceCas9Parameters(yamlFile)
#' is_casDataRun()
#' # this is a header to a document - not a live Snakemake workspace - demo
#' # how to actually run ad-hoc the code
#' bamFile <- system.file("extdata", "cas9_FAK76554.bam", package = "nanopoRe")
#' referenceGenome <- system.file("extdata", "cas9_demo_ref.fasta", package = "nanopoRe")
#' bedTargets <- system.file("extdata", "cas9_demo_target.bed", package = "nanopoRe")
#' unmappedQ <- system.file("extdata", "cas9_FAK76554.unmapped.quals", package = "nanopoRe")
#' setCas9ParameterValue("reference_genome", referenceGenome)
#' setCas9ParameterValue("target_regions", bedTargets)
#' addCas9ParameterValue("bam_file", bamFile)
#' addCas9ParameterValue("unmapped_quals", unmappedQ)
#' # lines between the previous comment and here are typically not required
#' RunEnrichmentAnalysis(mc.cores=1)
#'
#' @export
RunEnrichmentAnalysis <- function(force=FALSE, ...) {

    if (!is_casDataRun() || force) {
        message(paste0("loading reference genome: ",
            getCas9ParameterValue("reference_genome")))
        setReferenceGenome(getCas9ParameterValue("reference_genome"))
        loadReferenceGenome()

        # define BAM location
        if (!hasCas9ParameterField("bam_file")) {
            addCas9ParameterValue("bam_file", file.path("Analysis", "Minimap2",
                paste0(getCas9ParameterValue("study_name"), ".bam")))
        }
        if (!hasCas9ParameterField("unmapped_quals")) {
            addCas9ParameterValue("unmapped_quals", file.path("Analysis", "Minimap2", paste0(getCas9ParameterValue("study_name"), ".unmapped.quals")))
        }

        message(paste0("loading unmapped_reads: ", getCas9ParameterValue("unmapped_quals")))
        harvestUnmappedQuals(getCas9ParameterValue("unmapped_quals"), force)

        message(paste0("parsing genome coordinates: ", "..."))
        parseGenomeCoordinates(getCas9ParameterValue("target_regions"),
            getCas9ParameterValue("target_proximity"))

        message(paste0("parsing BAM file: ", getCas9ParameterValue("bam_file")))
        parseEnrichmentBam(
            getCas9ParameterValue("bam_file"),
            getCas9ParameterValue("gstride"),
            getCas9ParameterValue("offtarget_level"))
        message(paste0("preparing mapping characteristics", "\n"))
        mapUniverseTypes(...)

        dir.create(file.path(getRpath(), "..", "OnTarget"), showWarnings = FALSE, recursive=TRUE)
        aggregatedOntarget(...)

        dir.create(file.path(getRpath(), "..", "OffTarget"), showWarnings = FALSE, recursive=TRUE)
        aggregatedOfftarget(...)
    }

}


aggregatedOntarget <- function(mc.cores=min(parallel::detectCores()-1, max_threads)) {
    message(paste0("scoring data per target region", "\n"))
    ontargetUniverse <- get("ontargetUniverse", envir = getCachedObject("GRanges"))
    max_threads <- getCas9ParameterValue("threads")
    suppressWarnings({
        aggregatedCov <- dplyr::bind_rows(pbmclapply(seq_along(ontargetUniverse), aggregateDepthInfo, xr=ontargetUniverse, ontarget=TRUE, mc.cores=mc.cores), .id = "column_label")
        aggregatedCovFile <- file.path(getRpath(), paste0(getCas9ParameterValue("study_name"), "_aggregated_coverage", ".Rdata"))
        aggregatedGR <- GenomicRanges::makeGRangesFromDataFrame(aggregatedCov[,-1], keep.extra.columns = TRUE)
        save(aggregatedGR, file=aggregatedCovFile)
    })
}


aggregatedOfftarget <- function(mc.cores=min(parallel::detectCores()-1, max_threads)) {
    message(paste0("parsing off-target/background coverage - please be patient ...", "\n"))
    offtargetUniverse <- get("offtargetUniverse", envir = getCachedObject("GRanges"))
    max_threads <- getCas9ParameterValue("threads")
    suppressWarnings({
        offtCov <- pbmclapply(seq_along(offtargetUniverse), aggregateDepthInfo, xr=offtargetUniverse, geneId="OffTarget", mc.cores=mc.cores)
        aggregatedOff <- dplyr::bind_rows(offtCov, .id = "column_label")
        aggregatedOffFile <- file.path(getRpath(), paste0(getCas9ParameterValue("study_name"), "_aggregated_offt_coverage", ".Rdata"))
        save(aggregatedOff, file=aggregatedOffFile)
    })
}



#' @importFrom IRanges start<-
#' @importFrom IRanges end<-
#' @importFrom stats end
#' @importFrom GenomicRanges strand
aggregateDepthInfo <- function(x, xr, ontarget=FALSE, geneId=NA) {
    targetRegion <- xr[x]
    proximalRegion <- targetRegion
    target_proximity <- getCas9ParameterValue("target_proximity")
    if (ontarget) {
        start(proximalRegion) <- max(IRanges::start(proximalRegion) - target_proximity, 1)
        end(proximalRegion) <- IRanges::end(proximalRegion) + target_proximity
    }
    if (is.na(geneId)) {
        geneId <- names(targetRegion)
    }
    wga <- get("wga", envir = getCachedObject("GRanges"))
    referenceGenome <- get("referenceGenome", envir = get(getEnvironment()))
    referenceGenomeSequence <- get("referenceGenomeSequence", envir = get(getEnvironment()))
    c0 <- wga[seqnames(wga)==as.character(seqnames(targetRegion))]
    c1 <- c0[S4Vectors::subjectHits(GenomicRanges::findOverlaps(proximalRegion, GenomicRanges::granges(c0)))]
    c2 <- c0[S4Vectors::subjectHits(GenomicRanges::findOverlaps(targetRegion, GenomicRanges::granges(c0)))]
    seqlevels(c1) <- unique(as.character(seqnames(c1)))
    seqnames(c1) <- factor(seqnames(c1))
    cov <- GenomicAlignments::coverage(c1, shift=-IRanges::start(proximalRegion), width=IRanges::width(proximalRegion))
    seqlen <- IRanges::width(proximalRegion)
    names(seqlen)<-as.character(seqnames(targetRegion))
    bins <- GenomicRanges::tileGenome(seqlengths=seqlen, tilewidth=10, cut.last.tile.in.chrom=TRUE)
    ba <- binnedAverage(bins, cov, "binned_cov")
    slen <- IRanges::width(referenceGenomeSequence[getStringSetId(unique(seqnames(ba)))])
    names(slen) <- unique(seqnames(ba))
    seqlengths(ba) <- slen
    ba$pos <- IRanges::start(ba)
    IRanges::end(ba) <- IRanges::end(ba)+IRanges::start(targetRegion)
    IRanges::start(ba) <- IRanges::start(ba)+IRanges::start(targetRegion)
    ba$gene <- geneId

    # calculate the coverage of reads on fwd strand for directionality plotting
    fcov <- GenomicAlignments::coverage(c1[which(as.character(strand(c1))=="+")], shift=-IRanges::start(proximalRegion), width=IRanges::width(proximalRegion))
    ba$fwd_cov <- mcols(binnedAverage(bins, fcov, "fwd_cov"))$fwd_cov
    # write the target sequences to file ...
    if (ontarget) {
        write(mcols(c2)$qname, file.path(file.path(getRpath(), "..", "OnTarget"), paste0(getCas9ParameterValue("study_name"), ".", geneId, ".mappedreads")), append=TRUE)
    } else {
        write(mcols(c2)$qname, file.path(file.path(getRpath(), "..", "OffTarget"), paste0(getCas9ParameterValue("study_name"), ".", geneId, ".mappedreads")), append=TRUE)
    }
    return(as.data.frame(ba))
}






mapUniverseTypes <- function(...) {
    message(paste0("background", "\n"))
    backgroundUniverse <- bamMineUniverse(get("backgroundUniverse", envir = getCachedObject("GRanges")), ...)
    message(paste0("offtarget", "\n"))
    offtargetUniverse <- bamMineUniverse(get("offtargetUniverse", envir = getCachedObject("GRanges")), ...)
    message(paste0("ontarget", "\n"))
    ontargetUniverse <-bamMineUniverse(get("ontargetUniverse", envir = getCachedObject("GRanges")), ...)
    message(paste0("target proximal", "\n"))
    targetproximalUniverse <-bamMineUniverse(get("targetproximalUniverse", envir = getCachedObject("GRanges")), ...)

    br <- get("br", envir = getCachedObject("GRanges"))
    wga.cov <- get("wga.cov", envir = getCachedObject("GRanges"))

    mappingResultsFile <- file.path(getRpath(), paste0(getCas9ParameterValue("study_name"), "_mapping_results", ".Rdata"))
    message(paste0("writing result to:", mappingResultsFile))
    save(br, wga.cov, backgroundUniverse, offtargetUniverse, ontargetUniverse, targetproximalUniverse, file=mappingResultsFile)

    assign("backgroundUniverse", backgroundUniverse, envir = getCachedObject("GRanges"))
    assign("offtargetUniverse", offtargetUniverse, envir = getCachedObject("GRanges"))
    assign("targetproximalUniverse", targetproximalUniverse, envir = getCachedObject("GRanges"))
    assign("ontargetUniverse", ontargetUniverse, envir = getCachedObject("GRanges"))
}





parseEnrichmentBam <- function(mappedBam, gstride, offtarget_level) {
    if (!(exists("referenceGenome", envir = get(getEnvironment())) &&
        exists("referenceGenomeSequence", envir = get(getEnvironment())))) {
            loadReferenceGenome()
    }
    referenceGenome <- get("referenceGenome", envir = get(getEnvironment()))
    referenceGenomeSequence <- get("referenceGenomeSequence", envir = get(getEnvironment()))

    # set GRanges sequence lengths
    gr <- get("gr", envir = getCachedObject("GRanges"))
    br <- get("br", envir = getCachedObject("GRanges"))
    fr <- get("fr", envir = getCachedObject("GRanges"))
    seqlengths(gr) <- width(
        referenceGenomeSequence[getStringSetId(names(seqlengths(gr)))])

    message(paste0("primary mappings"))
    params=Rsamtools::ScanBamParam(
        which=gr, what=c("mapq", "qual", "qname"), tag=c("NM"),
        flag=Rsamtools::scanBamFlag(
            isSupplementaryAlignment=FALSE, isSecondaryAlignment=FALSE))
    wga <- GenomicAlignments::readGAlignments(file=mappedBam, param=params)

    max.read.len <- max(width(wga))
    wgac <- coverage(wga)
    bins <- tileGenome(seqlengths(gr), tilewidth=gstride, cut.last.tile.in.chrom=TRUE)
    # handle any hanger-on chromosome ids
    wgac <- wgac[which(names(wgac) %in% seqlevels(bins))]
    wgaba.all <- binnedAverage(bins, wgac, "binned_cov")

    message(paste0("looking for off-target"))
    wgaba <- wgaba.all[-S4Vectors::queryHits(GenomicRanges::findOverlaps(wgaba.all, br))]
    wgaba <- wgaba[-S4Vectors::queryHits(GenomicRanges::findOverlaps(wgaba, fr))]
    wga.cov <- mean(wgaba$binned_cov) # calculate mean background coverage
    offR <- wgaba[which(wgaba$binned_cov > (offtarget_level*wga.cov))]
    backgroundR <- wgaba[-S4Vectors::queryHits(GenomicRanges::findOverlaps(wgaba, offR))]

    assign("wga", wga, envir = getCachedObject("GRanges"))
    assign("backgroundUniverse", GenomicRanges::reduce(backgroundR), envir = getCachedObject("GRanges"))
    assign("offtargetUniverse", GenomicRanges::reduce(offR), envir = getCachedObject("GRanges"))
    assign("targetproximalUniverse", GenomicRanges::reduce(fr), envir = getCachedObject("GRanges"))
    ontargetUniverse <- GRanges(br)
    names(ontargetUniverse) <- names(br)
    assign("ontargetUniverse", ontargetUniverse, envir = getCachedObject("GRanges"))
    assign("wga.cov", wga.cov, envir = getCachedObject("GRanges"))
}



#' @importFrom stats end
#' @importFrom XVector subseq
getStartStrand <- function(x, gdata) {
    wga <- get("wga", envir = getCachedObject("GRanges"))
    referenceGenome <- get("referenceGenome", envir = get(getEnvironment()))
    referenceGenomeSequence <- get("referenceGenomeSequence", envir = get(getEnvironment()))
    nd <- gdata[x]
    # how to perform a fast reduction on data size - filter by chromosome
    c0 <- wga[seqnames(wga)==as.character(seqnames(nd))]
    c1 <- c0[S4Vectors::subjectHits(GenomicRanges::findOverlaps(nd, GenomicRanges::granges(c0)))]

    starts <- which(start(c1)>=start(nd))
    c2 <- c1[starts]

    seqlevels(c1) <- unique(as.character(seqnames(c1)))
    seqnames(c1) <- factor(seqnames(c1))

    cov <- GenomicAlignments::coverage(c1, shift=-start(nd), width=width(nd))
    depths <- rep(as.integer(unlist(S4Vectors::runValue(cov))), as.integer(unlist(S4Vectors::runLength(cov))))

    qdata <- quantile(depths)
    seqnm <- as.character(seqnames(nd))
    seqnmssid <- getStringSetId(seqnm)
    seqnmseq <- referenceGenomeSequence[[seqnmssid]]

    bss <- XVector::subseq(seqnmseq, start=IRanges::start(nd), end=IRanges::end(nd))
    lf <- Biostrings::letterFrequency(bss, c("G", "C", "N"))
    rstart <- length(c2)                                       # this is equiv. to earlier readStarts
    basesstart <- sum(GenomicAlignments::qwidth(c2))                              # earlier basesReadsStarted
    meanreadlen <- mean(GenomicAlignments::qwidth(c1))                            # see previous *readLen* - now mean for *any* overlapping read
    startreadlen <- mean(GenomicAlignments::qwidth(c2))                           # what is the mean sequence length for reads that start in segment?
    strandp <- length(which(as.character(GenomicAlignments::strand(c1))=="+"))                  # within the interval; not linked to readStarts
    strandn <- length(which(as.character(GenomicAlignments::strand(c1))=="-"))
    gccount <- lf[1]+lf[2]
    ncount <- lf[3]
    # discussion with Olle on how Qvalues are best summarised and prepared ... ShortRead alphabetScore is not ideal
    # while alphabetScore(mcols(c2)$qual) / width(mcols(c2)$qual)) gives a number this is not scaled appropriately
    mapq <- phredmean(mcols(c2)$mapq)
    map0 <- phredmean(mcols(c1)$mapq)
    readq <- phredmean(unlist(lapply(as.character(mcols(c2)$qual), qualToMeanQ)))
    read0 <- phredmean(unlist(lapply(as.character(mcols(c1)$qual), qualToMeanQ)))

    d005 <- qdata[1]
    d025 <- qdata[2]
    d050 <- qdata[3]
    d075 <- qdata[4]
    d095 <- qdata[5]
    dmean <- mean(depths)
    nm <- sum(mcols(c2)$NM)                                    # this is the #NM mismatch count; reads starting in segment
    cigardel <- sum(sum(width(GenomicAlignments::cigarRangesAlongPairwiseSpace(GenomicAlignments::cigar(c2), ops=c("D")))))
    cigarins <- sum(sum(width(GenomicAlignments::cigarRangesAlongPairwiseSpace(GenomicAlignments::cigar(c2), ops=c("I")))))
    cigarmapped <- sum(GenomicAlignments::cigarWidthAlongQuerySpace(GenomicAlignments::cigar(c2), after.soft.clipping=TRUE))
    return(c(rstart, basesstart, meanreadlen, startreadlen, strandp, strandn, gccount, ncount, mapq, map0, readq, read0, d005, d025, d050, d075, d095, dmean, nm, cigardel, cigarins, cigarmapped))
}


bamMineUniverse <- function(universe, mc.cores=min(parallel::detectCores()-1, max_threads)) {

    max_threads <- getCas9ParameterValue("threads")

    startStrand <-matrix(unlist(pbmclapply(seq_along(seqnames(universe)), getStartStrand, gdata=universe, mc.cores=mc.cores)), ncol=22, byrow=TRUE)
    universe$rstart <- startStrand[,1]
    universe$basesstart <- startStrand[,2]
    universe$meanreadlen <- startStrand[,3]
    universe$startreadlen <- startStrand[,4]
    universe$strandp <- startStrand[,5]
    universe$strandn <- startStrand[,6]
    universe$gccount <- startStrand[,7]
    universe$ncount <- startStrand[,8]
    universe$mapq <- startStrand[,9]
    universe$map0 <- startStrand[,10]
    universe$readq <- startStrand[,11]
    universe$read0 <- startStrand[,12]
    universe$d005 <- startStrand[,13]
    universe$d025 <- startStrand[,14]
    universe$d050 <- startStrand[,15]
    universe$d075 <- startStrand[,16]
    universe$d095 <- startStrand[,17]
    universe$dmean <- startStrand[,18]
    universe$nm <- startStrand[,19]
    universe$cigardel <- startStrand[,20]
    universe$cigarins <- startStrand[,21]
    universe$cigarmapped <- startStrand[,22]
    return(universe)
}




#' @importFrom GenomicRanges flank
#' @importFrom GenomicRanges union
#' @importFrom GenomicRanges disjoin
#' @importFrom GenomicRanges reduce
parseGenomeCoordinates <- function(bed_src, target_proximity) {
    bed <- data.table::fread(file=bed_src)
    # create a genomic ranges object for the bed file elements
    br <- GRanges(seqnames=unlist(bed[,1]), IRanges(start=unlist(bed[,2]),
        end=unlist(bed[,3])))
    # define flanking regions for the target-proximal analysis
    fr <- GenomicRanges::union(
        GenomicRanges::flank(br, width=target_proximity, start=TRUE),
        GenomicRanges::flank(br, width=target_proximity, start=FALSE))
    # create a genomic ranges object for the chromosomes
    gr <- getReferenceGenomeGR()
    # and define an off-target ranges ..
    or <- GenomicRanges::setdiff(
        GenomicRanges::disjoin(
            c(gr, GenomicRanges::reduce(GenomicRanges::union(br, fr)))),
                GenomicRanges::reduce(GenomicRanges::union(br, fr)))
    # name genes associated with BED annotations for subsequent display ...
    names(br) <- unlist(bed[,4])

    setCachedObject("GRanges", new.env())
    assign("bed", bed, envir = getCachedObject("GRanges"))
    assign("br", br, envir = getCachedObject("GRanges"))
    assign("fr", fr, envir = getCachedObject("GRanges"))
    assign("gr", gr, envir = getCachedObject("GRanges"))
    assign("or", or, envir = getCachedObject("GRanges"))
}





harvestUnmappedQuals <- function(qualfilelocation, force=FALSE, chunk.size=100000) {

    chromosomeFile <- file.path(getRpath(), paste(sub("\\.[^.]*$", "", basename(qualfilelocation)), "rcounts", "Rdata",sep="."))
    if (file.exists(chromosomeFile) && !force) {
        unmapped.content <- readRDS(file=chromosomeFile)
        return(invisible(unmapped.content))
    }
    offset <- 0
    unmapped.content <- data.frame(width=integer(),
        quality=numeric(), stringsAsFactors=FALSE)
    repeat {
        message(paste("iter", offset))
        qual.data <- data.table::fread(file=qualfilelocation, nrows=chunk.size, skip=offset, sep="\n", header=FALSE)
        colnames(qual.data) <- "qual"
        width <- nchar(qual.data$qual)
        meanQ <- unlist(lapply(qual.data$qual, qualToMeanQ))
        unmapped.content <- rbind(unmapped.content, data.frame(width=width, quality=meanQ))
        if (nrow(qual.data) < chunk.size) {
            break
        }
        offset <- offset + chunk.size
    }
    saveRDS(unmapped.content, file=chromosomeFile)
    return(invisible(unmapped.content))

}
