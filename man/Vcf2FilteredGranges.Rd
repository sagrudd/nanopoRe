% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SnifflesSV.R
\name{Vcf2FilteredGranges}
\alias{Vcf2FilteredGranges}
\title{extract SV variants from a sniffles VCF file by SVTYPE}
\usage{
Vcf2FilteredGranges(vcfFile, svtype = "INS")
}
\arguments{
\item{vcfFile}{is the path to the vcfFile to plot}

\item{svtype}{defines one of INS/DEL/DUP}
}
\value{
a GRanges object describing the duplications
}
\description{
This method will extract SV variants as called by Sniffles
into a GRanges object with the variation length declared in the
SVLEN field, read support in RE, support type in SUPTYPE. The
variants will be filtered according to the svtype parameter
}
\examples{
vcfFile <- system.file("extdata", "GM24385.nf7.chr20.vcf",
    package = "nanopoRe")
Vcf2FilteredGranges(vcfFile, svtype = "INS")

}
