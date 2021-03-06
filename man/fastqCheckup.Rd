% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FileHandling.R
\name{fastqCheckup}
\alias{fastqCheckup}
\title{perform a sanity check and prepare summary info on fastq file}
\usage{
fastqCheckup(input_fastq, force = FALSE)
}
\arguments{
\item{input_fastq}{is the path to the fastq file to check}

\item{force}{logical value describing whether the analysis should be force
recalculated}
}
\value{
vector of observations
}
\description{
This method will prepare summary statistics on a fastq file and cache
results (assuming immutable nature)
}
\examples{
fastqFile <- system.file("extdata", "example.fastq.gz", package = "nanopoRe")
fqcheck <- fastqCheckup(fastqFile)
names(fqcheck)
fqcheck[['reads']]

}
