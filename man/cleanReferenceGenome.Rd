% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReferenceGenome.R
\name{cleanReferenceGenome}
\alias{cleanReferenceGenome}
\title{cleanup the reference genome data}
\usage{
cleanReferenceGenome(delIds)
}
\arguments{
\item{delIds}{is a vector of chromosomes to be dropped}
}
\value{
None
}
\description{
This method will remove unwanted chromosomes from the reference genome
collection; useful for when analyses are to
be linked to a set of core autosomes, for example.
}
\examples{
init()
referenceFasta <- system.file("extdata",
    "Escherichia_coli_complete_genome.fasta",
    package = "nanopoRe")
setReferenceGenome(referenceFasta)
getChromosomeIds()
cleanReferenceGenome(getChromosomeIds()[1])
getChromosomeIds()

}
