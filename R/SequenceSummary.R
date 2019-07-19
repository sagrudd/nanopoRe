
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

  return(seqsumdata)
}

