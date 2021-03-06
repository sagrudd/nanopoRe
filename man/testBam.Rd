% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BamParser.R
\name{testBam}
\alias{testBam}
\title{extract content from a BAM file}
\usage{
testBam(bamFile, yieldSize = 100L)
}
\arguments{
\item{bamFile}{is the location to the BAM file to parse}

\item{yieldSize}{is the number of mapping observations to extract}
}
\value{
list of summary observations
}
\description{
This method will extract a block of BAM content from a given BAM file
}
\examples{
# define the path to a BAM file
demoBam <- system.file("extdata",
    "Ecoli_zymo_R10_filt_subs.bam",
    package = "nanopoRe")
bamdata <- testBam(demoBam, 10)

}
