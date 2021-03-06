% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnrichmentPresentation.R
\name{isEnrichmentAnalysisComplete}
\alias{isEnrichmentAnalysisComplete}
\title{accessory method to test whether the cas9 data indexing
has been performed}
\usage{
isEnrichmentAnalysisComplete(resultsDir = NULL)
}
\arguments{
\item{resultsDir}{the location of the folder where files should be pulled
from - this should be determined automatically from config.yaml}
}
\value{
TRUE or FALSE depending
}
\description{
The cas9 workflow involves a deep dive into the experimental data and
produces a couple of Rdata files that contain the experimental context and
information relating to on-target, off-target and other effects. This method
tests whether the pre-process analysis has been performed
}
\examples{
yamlFile <- system.file("extdata", "cas9_demo.yaml", package = "nanopoRe")
importConfigYAML(yamlFile=yamlFile)
# some demo data is stored with the package - import the data by supply the
# non-standard path
isEnrichmentAnalysisComplete(system.file("extdata", package="nanopoRe"))

}
