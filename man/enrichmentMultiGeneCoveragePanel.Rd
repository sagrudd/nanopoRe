% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnrichmentPresentation.R
\name{enrichmentMultiGeneCoveragePanel}
\alias{enrichmentMultiGeneCoveragePanel}
\title{prepare a plot of depth-of-coverage across each of the cas9-enrichment
targets included in study}
\usage{
enrichmentMultiGeneCoveragePanel(colMax = 4)
}
\arguments{
\item{colMax}{defines the maximum number of columns to use within the tiled
matrix}
}
\value{
ggplot2 format graph
}
\description{
prepares a plot of depth-of-coverage across each of the cas9-enrichment
targets included in study. This is rendered as a multi-tile plot; simple for
quick visualisation of the target regions
}
\examples{
yamlFile <- system.file("extdata", "cas9_demo.yaml", package = "nanopoRe")
importConfigYAML(yamlFile=yamlFile)
importEnrichmentResults(system.file("extdata", package="nanopoRe"))
enrichmentMultiGeneCoveragePanel()

}
