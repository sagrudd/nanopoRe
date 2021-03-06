% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cigar.R
\name{cigarR2Q}
\alias{cigarR2Q}
\title{maps reference coordinates back to the query sequence}
\usage{
cigarR2Q(poi, cigar)
}
\arguments{
\item{poi}{is a vector of numeric positions of interest}

\item{cigar}{is the CIGAR string to parse}
}
\value{
matrix of coordinates
}
\description{
This method will prepare a vector of query coordinates from
the supplied reference coordinates
}
\examples{
cigarStr <-
    "46S14M5D14M2D44M3I6M2I14M3I15M1D27M1D11M1I20M1I8M2I2M3I10S"
cigarR2Q(c(20,40,60,80,100), cigarStr)

}
