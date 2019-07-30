##########################################################################
### This is a collection of specific methods to an ONT_tutorial
##########################################################################



#' ont_tutorial_sv = summarise the mapping characteristics from a given BAM file
#'
#' this method prepares a kable format table summarising the key mapping characteritics from a provided
#' BAM file summary (from bamSummarise). This will describe the mapped sequence reads that are defined
#' as primary, secondary and supplementary and will overlay additional data that releases to read
#' length, read quality and other characteristics as reported in the tutorial document
#'
#' @importFrom knitr kable
#' @importFrom tibble add_row
#' @importFrom kableExtra footnote
#' @importFrom kableExtra footnote_marker_symbol
#' @param bamSummary is the result from  bamSummarise
#' @param validationResponse is the product of the floundeR fastqCheckup method
#' @param bamFile  is a path to a file
#' @return kable table of mapping characteristics as tuned for the ont_tutorial_sv document
#'
#' @examples
#' \dontrun{
#' mappingFile <- system.file("extdata", "mapping.bam", package = "nanopoRe", mustWork = TRUE)
#' bamSummary <- bamSummarise(bamFile)
#' kable <- SVMappingCharacteristicTable(bamSummary)
#' }
#'
#' @export
SVMappingCharacteristicTable <- function(bamSummary, validationResponse, bamFile) {

  primaryBAM <- bamSummary[which(bamSummary$readFlag=="Primary"),]
  unmappedBAM <- bamSummary[which(bamSummary$readFlag=="Unmapped"),]

  totalSeqs <- as.numeric(validationResponse[7])
  primaryMaps <- nrow(primaryBAM)
  unmapped <- nrow(unmappedBAM)

  totalNucleotides <- as.numeric(validationResponse[6])
  mappedReadNt <- sum(primaryBAM$qwidth)
  clippedReadNt <- sum(primaryBAM$coverage * primaryBAM$qwidth)  # clipped_width ...
  unmappedNt <- sum(unmappedBAM$qwidth)

  mapping <- data.frame(key=character(0), value=character(0), percentage=character(0))

  # statistics on number of reads
  mapping <- mapping %>% add_row(key="fastq reads", value=scales::comma_format()(totalSeqs), percentage="")

  mapping <- mapping %>% add_row(key="mapped reads (primary)",
                                 value=scales::comma_format()(primaryMaps),
                                 percentage=paste0(round(primaryMaps / totalSeqs * 100, digits=2),"%"))
  mapping <- mapping %>% add_row(key="mapped reads (secondary)",
                                 value=scales::comma_format()(length(which(bamSummary$readFlag=="Secondary"))),
                                 percentage="")
  mapping <- mapping %>% add_row(key="mapped reads (supplementary)",
                                 value=scales::comma_format()(length(which(bamSummary$readFlag=="Supplementary"))),
                                 percentage="")
  mapping <- mapping %>% add_row(key="unmapped reads",
                                 value=scales::comma_format()(unmapped),
                                 percentage=paste0(round(unmapped / totalSeqs * 100, digits=2),"%"))

  mapping <- mapping %>% add_row(key="fastq nucleotides", value=scales::comma_format()(totalNucleotides), percentage="")

  mapping <- mapping %>% add_row(key="nucleotides from primary mapped reads",
                                 value=scales::comma_format()(mappedReadNt),
                                 percentage=paste0(round(mappedReadNt / totalNucleotides * 100, digits=2),"%"))
  mapping <- mapping %>% add_row(key="clipped nucleotides from primary mappings",
                                 value=scales::comma_format()(clippedReadNt),
                                 percentage=paste0(round(clippedReadNt / totalNucleotides * 100, digits=2),"%"))
  mapping <- mapping %>% add_row(key="unmapped nucleotides",
                                 value=scales::comma_format()(unmappedNt),
                                 percentage=paste0(round(unmappedNt / totalNucleotides * 100, digits=2),"%"))
  mapping <- mapping %>% add_row(key="nucleotides mapped (secondary mappings)",
                                 value=scales::comma_format()(sum(bamSummary[which(bamSummary$readFlag=="Secondary"),]$qwidth * bamSummary[which(bamSummary$readFlag=="Secondary"),]$coverage)),
                                 percentage="")
  mapping <- mapping %>% add_row(key="nucleotides mapped (supplementary mappings)",
                                 value=scales::comma_format()(sum(bamSummary[which(bamSummary$readFlag=="Supplementary"),]$qwidth * bamSummary[which(bamSummary$readFlag=="Supplementary"),]$coverage)),
                                 percentage="")

  mapping <- mapping %>% add_row(key="mean read length (primary mapping)",
                                 value=scales::comma_format()(mean(bamSummary[which(bamSummary$readFlag=="Primary"),]$qwidth)), percentage="")
  mapping <- mapping %>% add_row(key="mean read length (secondary mapping)",
                                 value=scales::comma_format()(mean(bamSummary[which(bamSummary$readFlag=="Secondary"),]$qwidth)), percentage="")
  mapping <- mapping %>% add_row(key="mean read length (supplementary mapping)",
                                 value=scales::comma_format()(mean(bamSummary[which(bamSummary$readFlag=="Supplementary"),]$qwidth)), percentage="")
  mapping <- mapping %>% add_row(key="mean read length (unmapped)",
                                 value=scales::comma_format()(mean(bamSummary[which(bamSummary$readFlag=="Unmapped"),]$qwidth)), percentage="")


  mapping <- mapping %>% add_row(key="mean read quality (primary mapping)",
                                 value=paste0(round(phredmean(bamSummary[which(bamSummary$readFlag=="Primary"),]$readq), digits=2)), percentage="")
  mapping <- mapping %>% add_row(key="mean read quality (secondary mapping)",
                                 value=paste0(round(phredmean(bamSummary[which(bamSummary$readFlag=="Secondary"),]$readq), digits=2)), percentage="")
  mapping <- mapping %>% add_row(key="mean read quality (supplementary mapping)",
                                 value=paste0(round(phredmean(bamSummary[which(bamSummary$readFlag=="Supplementary"),]$readq), digits=2)), percentage="")
  mapping <- mapping %>% add_row(key="mean read quality (unmapped)",
                                 value=paste0(round(phredmean(bamSummary[which(bamSummary$readFlag=="Unmapped"),]$readq), digits=2)), percentage="")


  mapping <- mapping %>% add_row(key="mean mapping quality (primary mapping)",
                                 value=paste0(round(mean(bamSummary[which(bamSummary$readFlag=="Primary"),]$mapq), digits=2)), percentage="")
  mapping <- mapping %>% add_row(key="mean mapping quality (secondary mapping)",
                                 value=paste0(round(mean(bamSummary[which(bamSummary$readFlag=="Secondary"),]$mapq), digits=2)), percentage="")
  mapping <- mapping %>% add_row(key="mean mapping quality (supplementary)",
                                 value=paste0(round(mean(bamSummary[which(bamSummary$readFlag=="Supplementary"),]$mapq), digits=2)), percentage="")


  mapping <- mapping %>% add_row(key="mean coverage (primary mapping)",
                                 value=paste0(round(mean(bamSummaryToCoverage(bamFile, flag="Primary")$binned_cov), digits=2)),
                                 percentage="")
  mapping <- mapping %>% add_row(key="mean coverage (secondary mapping)",
                                 value=paste0(round(mean(bamSummaryToCoverage(bamFile, flag="Secondary")$binned_cov), digits=2)),
                                 percentage="")
  mapping <- mapping %>% add_row(key="mean coverage (supplementary)",
                                 value=paste0(round(mean(bamSummaryToCoverage(bamFile, flag="Supplementary")$binned_cov), digits=2)),
                                 percentage="")

  rownames(mapping) <- mapping[,1]
  mapping <- mapping[,-1]

  row.names(mapping)[7]<- paste0(row.names(mapping)[7], footnote_marker_symbol(1, "html"))
  row.names(mapping)[8]<- paste0(row.names(mapping)[8], footnote_marker_symbol(2, "html"))

  ktable <- kable(mapping, format="html", col.names=rep(" ", ncol(mapping)), caption="Table summarising mapping characteristics", booktabs=TRUE, table.envir='table*', linesep="", escape = FALSE) %>%
    kable_styling(c("striped", "condensed"))  %>%
    pack_rows("Sequence Read Mapping", 1, 5) %>%
    pack_rows("Nucleotides Mapped", 6, 11) %>%
    pack_rows("Reads lengths", 12, 15) %>%
    pack_rows("Reads Qualities", 16, 19) %>%
    pack_rows("Mapping Qualities", 20, 22) %>%
    pack_rows("Depth of Coverage", 23, 25) %>%
    footnote(symbol=c("fastq bases are calculated from the qwidth field of the primary mapped sequences", "soft-clipped nucleotide bases are calculated from the BAM CIGAR annotation"), symbol_title="please note: ", footnote_as_chunk = TRUE)

  return(ktable)
}
