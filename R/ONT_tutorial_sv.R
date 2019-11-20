### This is a collection of specific methods to an ONT_tutorial



getBAMReadMappingStats <- function(totalSeqs, primaryMaps, bamSummary) {
    mapping <- data.frame(key = character(0), value = character(0),
        percentage = character(0))
    mapping <- mapping %>% add_row(key = "fastq reads", value = (
        scales::comma_format())(totalSeqs), percentage = "")
    mapping <- mapping %>% add_row(key = "mapped reads (primary)", value = (
        scales::comma_format())(primaryMaps), percentage = "")
    mapping <- mapping %>% add_row(key = "mapped reads (secondary)", value = (
        scales::comma_format())(length(which(bamSummary$readFlag ==
        "Secondary"))), percentage = "")
    mapping <- mapping %>% add_row(key="mapped reads (supplementary)", value=(
        scales::comma_format())(length(which(bamSummary$readFlag ==
        "Supplementary"))), percentage = "")
    return(mapping)
}



getBAMReadNucleotideStats <- function(totalNucleotides, mappedReadNt,
    clippedReadNt, bamSummary) {
    mapping <- data.frame(key = character(0), value = character(0),
        percentage = character(0))

    mapping <- mapping %>% add_row(key = "fastq nucleotides", value=(
        scales::comma_format())(totalNucleotides), percentage = "")
    mapping <- mapping %>% add_row(key="nucleotides from primary mapped reads",
        value = (scales::comma_format())(mappedReadNt), percentage = paste0(
        round(mappedReadNt/totalNucleotides * 100, digits = 2), "%"))
    mapping <- mapping %>% add_row(
        key="clipped nucleotides from primary mappings",
        value = (scales::comma_format())(clippedReadNt), percentage = paste0(
            round(clippedReadNt/totalNucleotides * 100, digits = 2), "%"))
    mapping <- mapping %>% add_row(
        key="nucleotides mapped (secondary mappings)", value = (
        scales::comma_format())(sum(bamSummary[which(bamSummary$readFlag ==
        "Secondary"), ]$qwidth * bamSummary[which(bamSummary$readFlag ==
        "Secondary"), ]$coverage)), percentage = "")
    mapping <- mapping %>% add_row(
        key = "nucleotides mapped (supplementary mappings)",
        value=(scales::comma_format())(sum(
        bamSummary[which(bamSummary$readFlag =="Supplementary"), ]$qwidth *
        bamSummary[which(bamSummary$readFlag == "Supplementary"), ]$coverage)),
        percentage = "")
    return(mapping)
}



getBAMReadLengthStats <- function(bamSummary) {
    mapping <- data.frame(
        key = character(0), value = character(0), percentage = character(0))

    mapping <- mapping %>% add_row(key = "mean read length (primary mapping)",
        value = (scales::comma_format())(mean(bamSummary[which(
        bamSummary$readFlag == "Primary"), ]$qwidth)), percentage = "")
    mapping <- mapping %>% add_row(key ="mean read length (secondary mapping)",
        value = (scales::comma_format())(mean(bamSummary[which(
        bamSummary$readFlag == "Secondary"), ]$qwidth)), percentage = "")
    mapping <- mapping %>% add_row(key=
        "mean read length (supplementary mapping)", value = (
        scales::comma_format())(mean(bamSummary[which(bamSummary$readFlag ==
        "Supplementary"), ]$qwidth)), percentage = "")

    return(mapping)
}



getBAMReadQualityStats <- function(bamSummary) {
    mapping <- data.frame(
        key = character(0), value = character(0), percentage = character(0))

    mapping <- mapping %>% add_row(
        key = "mean read quality (primary mapping)",
        value = paste0(round(phredmean(bamSummary[which(bamSummary$readFlag ==
        "Primary"), ]$readq), digits = 2)), percentage = "")
    mapping <- mapping %>% add_row(
        key = "mean read quality (secondary mapping)",
        value = paste0(round(phredmean(bamSummary[which(bamSummary$readFlag ==
        "Secondary"), ]$readq), digits = 2)), percentage = "")
    mapping <- mapping %>% add_row(
        key = "mean read quality (supplementary mapping)",
        value = paste0(round(phredmean(bamSummary[which(bamSummary$readFlag ==
        "Supplementary"), ]$readq), digits = 2)), percentage = "")
    mapping <- mapping %>% add_row(
        key = "mean mapping quality (primary mapping)",
        value = paste0(round(mean(bamSummary[which(bamSummary$readFlag ==
        "Primary"), ]$mapq), digits = 2)), percentage = "")
    mapping <- mapping %>% add_row(
        key = "mean mapping quality (secondary mapping)",
        value = paste0(round(mean(bamSummary[which(bamSummary$readFlag ==
        "Secondary"), ]$mapq), digits = 2)), percentage = "")
    mapping <- mapping %>% add_row(
        key = "mean mapping quality (supplementary)",
        value = paste0(round(mean(bamSummary[which(bamSummary$readFlag ==
        "Supplementary"), ]$mapq), digits = 2)), percentage = "")

    return(mapping)
}




getBAMReadCoverageStats <- function(bamFile, bamSummary) {
    mapping <- data.frame(
        key = character(0), value = character(0), percentage = character(0))

    mapping <- mapping %>% add_row(key = "mean coverage (primary mapping)",
        value = paste0(round(mean(bamSummaryToCoverage(bamFile,
        flag = "Primary")$binned_cov), digits = 2)), percentage = "")
    mapping <- mapping %>% add_row(key = "mean coverage (secondary mapping)",
        value = paste0(round(mean(bamSummaryToCoverage(bamFile,
        flag = "Secondary")$binned_cov), digits = 2)), percentage = "")
    mapping <- mapping %>% add_row(key = "mean coverage (supplementary)",
        value = paste0(round(mean(bamSummaryToCoverage(bamFile,
        flag = "Supplementary")$binned_cov), digits = 2)), percentage = "")
    return(mapping)
}

#' ont_tutorial_sv = summarise the mapping characteristics
#' from a given BAM file
#'
#' this method prepares a kable format table summarising the key mapping
#' characteritics from a provided BAM file summary (from bamSummarise). This
#' will describe the mapped sequence reads that are defined as primary,
#' secondary and supplementary and will overlay additional data that releases
#' to read length, read quality and other characteristics as reported in the
#' tutorial document
#'
#' @importFrom knitr kable
#' @importFrom tibble add_row
#' @importFrom kableExtra footnote
#' @importFrom kableExtra footnote_marker_symbol
#' @importFrom kableExtra group_rows pack_rows
#' @param bamSummary is the result from  bamSummarise
#' @param validationResponse is the product of the floundeR fastqCheckup method
#' @param bamFile  is a path to a file
#' @param format format of table to prepare
#' @return kable table of mapping characteristics as tuned for the
#' ont_tutorial_sv document
#'
#' @examples
#' demoBam <- system.file("extdata",
#'     "Ecoli_zymo_R10_filt_subs.bam",
#'     package = "nanopoRe")
#' demoFastq <- system.file("extdata",
#'     "Ecoli_zymo_R10.fastq.gz",
#'     package = "nanopoRe")
#' referenceGenome <- system.file("extdata",
#'     "Escherichia_coli_complete_genome.fasta",
#'     package = "nanopoRe")
#' setReferenceGenome(referenceGenome)
#' loadReferenceGenome()
#' bamSummary <- bamSummarise(demoBam, force=FALSE, blockSize=1000L)
#' fastqView <- fastqCheckup(demoFastq)
#' kable <- SVMappingCharacteristicTable(bamSummary, fastqView, demoBam)
#'
#' @export
SVMappingCharacteristicTable <- function(
    bamSummary, validationResponse, bamFile, format=NA) {
    primaryBAM <- bamSummary[which(bamSummary$readFlag == "Primary"), ]
    unmappedBAM <- bamSummary[which(bamSummary$readFlag == "Unmapped"), ]
    totalSeqs <- as.numeric(validationResponse[7])
    primaryMaps <- nrow(primaryBAM)
    unmapped <- nrow(unmappedBAM)
    totalNucleotides <- as.numeric(validationResponse[6])
    mappedReadNt <- sum(primaryBAM$qwidth)
    clippedReadNt <- sum(primaryBAM$coverage * primaryBAM$qwidth)
    unmappedNt <- sum(unmappedBAM$qwidth)

    mapping <- getBAMReadMappingStats(totalSeqs, primaryMaps, bamSummary)
    mapping <- rbind(mapping, getBAMReadNucleotideStats(totalNucleotides,
        mappedReadNt, clippedReadNt, bamSummary))
    mapping <- rbind(mapping, getBAMReadLengthStats(bamSummary))
    mapping <- rbind(mapping, getBAMReadQualityStats(bamSummary))
    mapping <- rbind(mapping, getBAMReadCoverageStats(bamFile, bamSummary))

    rownames(mapping) <- mapping[, 1]
    mapping <- mapping[, -1]

    if (is.na(format)) {
        return(mapping)
    } else if (format=="markdown") {
        simple <- kable(mapping, format = "markdown",
              caption = "Table summarising mapping characteristics",
              table.envir = "table*", linesep = "", escape = FALSE) %>%
            kable_styling(c("striped", "condensed"))
        return(simple)
    }

    row.names(mapping)[6] <- paste0(row.names(mapping)[6],
        footnote_marker_symbol(1, "html"))
    row.names(mapping)[7] <- paste0(row.names(mapping)[7],
        footnote_marker_symbol(2, "html"))
    ktable <- kable(mapping, format = "html", col.names=rep(" ",ncol(mapping)),
        caption = "Table summarising mapping characteristics", booktabs = TRUE,
        table.envir = "table*", linesep = "", escape = FALSE) %>%
        kable_styling(c("striped", "condensed")) %>%
        pack_rows("Read Mapping", 1, 4) %>%
        pack_rows("Nucleotide Mapping", 5, 9) %>%
        pack_rows("Read lengths", 10, 12) %>%
        pack_rows("Read Qualities", 13, 15) %>% pack_rows("Mapping Qualities",
        16, 18) %>% pack_rows("Depth of Coverage", 19, 21) %>%
        footnote(symbol = c(paste0("fastq bases are calculated from the qwidth",
        "field of the primary mapped sequences"),paste0("soft-clipped ",
        "nucleotide bases are calculated from the BAM CIGAR annotation")),
        symbol_title = "please note: ",
        footnote_as_chunk = TRUE)
    return(ktable)
}
