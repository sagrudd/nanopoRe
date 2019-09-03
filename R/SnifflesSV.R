

#' extract SV variants from a sniffles VCF file by SVTYPE
#'
#' This method will extract SV variants as called by Sniffles
#' into a GRanges object with the variation length declared in the
#' SVLEN field, read support in RE, support type in SUPTYPE. The
#' variants will be filtered according to the svtype parameter
#'
#' @importFrom vcfR read.vcfR
#' @importFrom vcfR extract.info
#' @importFrom vcfR getFIX
#' @param vcfFile is the path to the vcfFile to plot
#' @param svtype defines one of INS/DEL/DUP
#' @return a GRanges object describing the duplications
#'
#' @examples
#' \dontrun{
#' Vcf2FilteredGranges(vcfFile)
#' }
#'
#' @export
Vcf2FilteredGranges <- function(vcfFile, svtype = "INS") {
    # DEL DUP INS

    variants <- Vcf2GRanges(vcfFile)
    SVTYPE <- variants$SVTYPE
    keys <- which(SVTYPE == svtype)

    return(variants[keys, ])
}


#' convert VCF content into a GRanges object for nanopoRe usage
#'
#' This method will prepare a karyogram of annotated SVs against the genome reference
#'
#' @importFrom IRanges start
#' @importFrom IRanges ranges
#' @importFrom IRanges ranges<-
#' @importFrom vcfR read.vcfR
#' @importFrom vcfR extract.info
#' @importFrom vcfR getFIX
#' @param vcfFile is the path to the vcfFile to plot
#' @return a ggbio/ggplot2 graph
#'
#' @examples
#' \dontrun{
#' Vcf2GRanges(vcfFile)
#' }
#'
#' @export
Vcf2GRanges <- function(vcfFile) {

    vcf <- read.vcfR(vcfFile)
    karyo <- GRanges(seqnames = as.factor(getFIX(vcf)[, "CHROM"]), ranges = IRanges(start = as.numeric(getFIX(vcf)[,
        "POS"]), end = as.numeric(getFIX(vcf)[, "POS"])))

    karyo$SVTYPE <- factor(extract.info(vcf, element = "SVTYPE"))
    karyo$SVLEN <- as.numeric(extract.info(vcf, element = "SVLEN"))
    karyo$RE <- as.numeric(extract.info(vcf, element = "RE"))
    karyo$SU <- as.character(extract.info(vcf, element = "SUPTYPE"))

    # fix the negative values for deletions
    delets <- which(karyo$SVLEN < 0)
    karyo$SVLEN[delets] <- karyo$SVLEN[delets] * -1
    ranges(karyo)[delets] <- IRanges(start = start(karyo)[delets], end = start(karyo)[delets] + karyo$SVLEN[delets])

    return(karyo)
}


#' prepare karyogram of annotated SVs
#'
#' This method will prepare a karyogram of annotated SVs against the genome reference
#'
#' @importFrom ggbio autoplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 scale_fill_brewer
#' @importFrom ggplot2 scale_color_brewer
#' @param vcfFile is the path to the vcfFile to plot
#' @return a ggbio/ggplot2 graph
#'
#' @examples
#' \dontrun{
#' snifflesKaryogram(vcfFile)
#' }
#'
#' @export
snifflesKaryogram <- function(vcfFile) {

    karyo <- Vcf2GRanges(vcfFile)

    karyo <- karyo[grep("(MT)|(\\..+)", as.character(seqnames(karyo)), invert = TRUE), ]
    factoids <- gtools::mixedsort(unique(as.character(seqnames(karyo))))
    seqlevels(karyo) <- factoids
    seqlengths(karyo) <- getSeqLengths(names(seqlengths(karyo)))

    if (length(unique(as.character(seqnames(karyo)))) == 1) {
        # this is a tutorial workflow ... may require debugging depending on how script is used in real
        # workflows?
        if (unique(gtools:::mixedsort(as.character(seqnames(karyo)))) == "4") {
            seqlevels(karyo) <- append(seq(1, 22), c("X", "Y"))
            seqlengths(karyo) <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
                101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895,
                57227415)
        }
    }

    karyogram <- autoplot(karyo, layout = "karyogram", aes_string(color = "SVTYPE", fill = "SVTYPE"),
        main = "Karyogram showing location and type of structural variations") + scale_fill_brewer(palette = "Set1") +
        scale_color_brewer(palette = "Set1")

    return(karyogram)


}




#' prepare figure of SV length distribution
#'
#' This method will prepare figure of SV length distribution
#'
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 scale_fill_brewer
#' @importFrom ggplot2 scale_color_brewer
#' @param vcfFile is the path to the vcfFile to plot
#' @return a ggplot2 graph
#'
#' @examples
#' \dontrun{
#' svLengthDistribution(vcfFile)
#' }
#'
#' @export
svLengthDistribution <- function(vcfFile) {
    sdata <- Vcf2GRanges(vcfFile)
    plot <- ggplot(as.data.frame(mcols(sdata)), aes_string(x = "SVLEN", fill = "SVTYPE")) + geom_histogram(bins = 70) +
        scale_x_log10() + scale_fill_brewer(palette = "Set1") + labs(title = "Plot showing frequency of log10-scaled SV lengths") +
        xlab("log10 sequence length (bases)") + ylab("Frequency")
    return(plot)
}
