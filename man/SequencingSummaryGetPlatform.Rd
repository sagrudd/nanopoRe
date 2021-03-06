% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SequencingSummary.R
\name{SequencingSummaryGetPlatform}
\alias{SequencingSummaryGetPlatform}
\title{identify the most likely sequencing platform used to create
the summary data}
\usage{
SequencingSummaryGetPlatform(seqsum = NA)
}
\arguments{
\item{seqsum}{is the data.frame object as prepared by importSequencingSummary}
}
\value{
platform is the nanopore platform [MinION/Flongle/PromethION]
}
\description{
when provided with the seqsum object from a Sequencing_summary.txt file
identify the
most likely sequencing platform used
}
\examples{
init()
seqsumFile <- system.file("extdata",
    "sequencing_summary.txt.bz2", package = "nanopoRe")
seqsum <- importSequencingSummary(seqsumFile)
platform <- SequencingSummaryGetPlatform(seqsum)
platform

}
