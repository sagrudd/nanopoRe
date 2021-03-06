% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SequencingSummary.R
\name{SequencingSummaryReadLengthQualityDensity}
\alias{SequencingSummaryReadLengthQualityDensity}
\title{plot a density map of sequence lengths and quality scores}
\usage{
SequencingSummaryReadLengthQualityDensity(seqsum=NA,
    binFilter = 2,
    qcThreshold = 7
)
}
\arguments{
\item{seqsum}{is the data.frame object as prepared by importSequencingSummary}

\item{binFilter}{is the minimum number of reads to include a cell in plot
(removes speckle)}

\item{qcThreshold}{is the QC threshold used}
}
\value{
ggplot2 showing densities of read length and quality distribution
}
\description{
plots the density plot of read length against read mean quality scores
}
\examples{
init()
seqsumFile <- system.file("extdata",
    "sequencing_summary.txt.bz2", package = "nanopoRe")
seqsum <- importSequencingSummary(seqsumFile)
plot <- SequencingSummaryReadLengthQualityDensity(seqsum)

}
