% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SequenceCharacteristics.R
\name{ncalc}
\alias{ncalc}
\title{calculate Nstatistics (e.g. N50) for given sequence collection}
\usage{
ncalc(len.vector, n)
}
\arguments{
\item{len.vector}{is a vector of sequence lengths}

\item{n}{is a numeric of the value to use}
}
\value{
numeric of the corresponding N value
}
\description{
Given a set of contigs, the N50 is defined as the sequence length of the
shortest contig at 50% of the total genome length
(https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics)
}
\examples{
# to calculate e.g. N50
ncalc(seq(1,10), 0.5)

}
