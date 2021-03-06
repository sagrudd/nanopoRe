% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnrichmentPresentation.R
\name{enrichmentTargetPerformanceTable}
\alias{enrichmentTargetPerformanceTable}
\title{summarise cas9 target performance}
\usage{
enrichmentTargetPerformanceTable()
}
\value{
kable table in html format
}
\description{
prepare tabular data showing the enrichment performance for the cas9 targets
included within an enrichment study
}
\examples{
yamlFile <- system.file("extdata", "cas9_demo.yaml", package = "nanopoRe")
importConfigYAML(yamlFile=yamlFile)
importEnrichmentResults(system.file("extdata", package="nanopoRe"))
enrichmentTargetPerformanceTable()

}
