
#' calculates the fractional number of bases according to supplied b parameter
#'
#' an accessory method for various logicals; simple fractional base calculator
#'
#' @importFrom fastmatch fmatch
#' @importFrom R.utils bunzip2
#' @importFrom R.utils gunzip
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param barcodeFile pointer to a barcode file as produced by Guppy
#' @return a data.frame with barcode information
#'
#' @examples
#' seqsumFile <- system.file(
#'     'extdata', 'sequencing_summary.txt.bz2', package = 'nanopoRe')
#' importSequencingSummary(seqsumFile)
#' barcodeFile <- system.file(
#'     'extdata', 'barcoding_summary.txt.bz2', package = 'nanopoRe')
#' SequencingSummaryBarcodeMerge(barcodeFile=barcodeFile)
#'
#' @export
SequencingSummaryBarcodeMerge <- function(seqsum=NA, barcodeFile=NULL) {

    seqsum <- handleSeqSumCache(seqsum)
    seqsum <- seqsum[which(seqsum$passes_filtering), ]

    if (!"barcode_arrangement" %in% colnames(seqsum)) {

        if (!is.null(barcodeFile) && file.exists(barcodeFile)) {
            barcodedata <- data.table::fread(barcodeFile, select = c(
                "read_id", "barcode_arrangement"),
                showProgress = TRUE, stringsAsFactors = FALSE)
            pso <- order(seqsum$read_id, method = "radix")
            seqsum <- seqsum[pso, ]

            bco <- order(barcodedata$read_id, method = "radix")
            barcodedata <- barcodedata[bco, ]

            barcodeMapping <- fmatch(seqsum$read_id, barcodedata$read_id)
            seqsum$barcode_arrangement <- barcodedata[barcodeMapping, c("barcode_arrangement")]
        }
    }

    setCachedObject("barcodedata", seqsum)
    return(invisible(seqsum))
}


#' presents an emojifont based infographic for barcode characteristics
#'
#' an accessory method for various logicals; simple fractional base calculator
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param bcthreshold the threshold number of reads for a barcode to be
#' considered (150)
#' @return a numeric value expressed in gigabases
#'
#' @examples
#' seqsumFile <- system.file(
#'     'extdata', 'sequencing_summary.txt.bz2', package = 'nanopoRe')
#' importSequencingSummary(seqsumFile)
#' barcodeFile <- system.file(
#'     'extdata', 'barcoding_summary.txt.bz2', package = 'nanopoRe')
#' SequencingSummaryBarcodeMerge(barcodeFile=barcodeFile)
#' SequenceSummaryBarcodeInfoGraphic()
#'
#' @export
SequenceSummaryBarcodeInfoGraphic <- function(seqsum=NA, bcthreshold = 150) {

    seqsum <- handleBarcodedCache(seqsum)

    barcodes <- 0
    infographicFile <- NULL

    if ("barcode_arrangement" %in% names(seqsum)) {
        barcodedata = plyr::count(seqsum$barcode_arrangement)
        barcodedata = subset(barcodedata, barcodedata$freq > bcthreshold)
        names(barcodedata) <- gsub("x", "barcode", names(barcodedata))
        if ("unclassified" %in% barcodedata$barcode) {
            barcodes <- nrow(
                barcodedata[-which(barcodedata$barcode == "unclassified"), ])
            barcodeUnass <- sum(barcodedata[-which(barcodedata$barcode == "unclassified"), "freq"])/sum(barcodedata$freq) *
                100
            barcodeRange <- range(subset(barcodedata, "barcode" != "unclassified")$freq)
        } else {
            barcodes <- nrow(barcodedata)
            barcodeUnass = 100
            barcodeRange <- range(barcodedata$freq)
        }
    }

    if (barcodes > 0) {
        infographicFile <- infoGraphicPlot3(identifier = "seqsumBarcodes", panelA = c(key = "Reads with barcode (%)",
            value = round(barcodeUnass, digits = 1), icon = "fa-pie-chart"), panelB = c(key = "Barcodes identified",
            value = barcodes, icon = "fa-barcode"), panelC = c(key = "barcode variance", value = paste(barcodeRange,
            collapse = "\n"), icon = "fa-sliders"))
    }

    return(infographicFile)
}



barcodeSeqSummary <- function(barcodeId, myBarcode, myVector, myMethod, xlist = NA) {
    subVector <- myVector[which(myBarcode == barcodeId)]
    params <- list(subVector)
    if (!is.na(xlist)) {
        params <- append(params, xlist)
    }
    do.call(myMethod, params)
}


#' tabulate information on the SequencingSummary barcode fields and status
#'
#' tabulate information on the SequencingSummary barcode fields and status
#'
#' @importFrom kableExtra kable_styling
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param bcthreshold the threshold number of reads for a barcode to be considered (150)
#' @return a kable table or NULL
#'
#' @examples
#' seqsumFile <- system.file('extdata', 'sequencing_summary.txt.bz2', package = 'nanopoRe')
#' importSequencingSummary(seqsumFile)
#' barcodeFile <- system.file('extdata', 'barcoding_summary.txt.bz2', package = 'nanopoRe')
#' SequencingSummaryBarcodeMerge(barcodeFile=barcodeFile)
#' SequenceSummaryBarcodeTable()
#'
#' @export
SequenceSummaryBarcodeTable <- function(seqsum=NA, bcthreshold = 150) {
    seqsum <- handleBarcodedCache(seqsum)

    barcodes <- 0
    barcodeTable <- NULL
    if ("barcode_arrangement" %in% names(seqsum)) {
        barcodedata = plyr::count(seqsum$barcode_arrangement)
        barcodedata = subset(barcodedata, barcodedata$freq > bcthreshold)
        names(barcodedata) <- gsub("x", "barcode", names(barcodedata))
        if ("unclassified" %in% barcodedata$barcode) {
            barcodes <- nrow(barcodedata[-which(barcodedata$barcode == "unclassified"), ])
            barcodeUnass <- sum(barcodedata[-which(barcodedata$barcode == "unclassified"), "freq"])/sum(barcodedata$freq) *
                100
            barcodeRange <- range(subset(barcodedata, "barcode" != "unclassified")$freq)
        } else {
            barcodes <- nrow(barcodedata)
            barcodeUnass = 100
            barcodeRange <- range(barcodedata$freq)
        }
    }
    if (barcodes > 0) {

        barcodedata <- cbind(barcodedata, `%` = round(barcodedata$freq/sum(barcodedata$freq) * 100,
            digits = 1))
        barcodedata <- cbind(barcodedata, Mb = round(unlist(lapply(as.character(barcodedata$barcode),
            barcodeSeqSummary, myBarcode = seqsum$barcode_arrangement, myVector = seqsum$sequence_length_template,
            myMethod = "sum"))/1e+06, digits = 0))
        barcodedata <- cbind(barcodedata, min = unlist(lapply(as.character(barcodedata$barcode), barcodeSeqSummary,
            myBarcode = seqsum$barcode_arrangement, myVector = seqsum$sequence_length_template, myMethod = "min")))
        barcodedata <- cbind(barcodedata, max = unlist(lapply(as.character(barcodedata$barcode), barcodeSeqSummary,
            myBarcode = seqsum$barcode_arrangement, myVector = seqsum$sequence_length_template, myMethod = "max")))
        barcodedata <- cbind(barcodedata, mean = round(unlist(lapply(as.character(barcodedata$barcode),
            barcodeSeqSummary, myBarcode = seqsum$barcode_arrangement, myVector = seqsum$sequence_length_template,
            myMethod = "mean")), digits = 0))
        barcodedata <- cbind(barcodedata, N50 = unlist(lapply(as.character(barcodedata$barcode), barcodeSeqSummary,
            myBarcode = seqsum$barcode_arrangement, myVector = seqsum$sequence_length_template, myMethod = "ncalc",
            xlist = list(n = 0.5))))
        barcodedata <- cbind(barcodedata, L50 = unlist(lapply(as.character(barcodedata$barcode), barcodeSeqSummary,
            myBarcode = seqsum$barcode_arrangement, myVector = seqsum$sequence_length_template, myMethod = "lcalc",
            xlist = list(n = 0.5))))
        barcodeTable <- kable(barcodedata, format = "html", caption = "Table summarising sequence collections ranked by barcode annotation",
            booktabs = TRUE, table.envir = "table*", linesep = "", escape = FALSE) %>% kable_styling(c("striped",
            "condensed"))
    }
    return(barcodeTable)
}





#' present a histogram of read count by sorted barcode Id
#'
#' present a ggplot2 histogram of read counts by sorted barcode Id
#'
#' @param seqsum is the data.frame object as prepared by importSequencingSummary
#' @param bcthreshold the threshold number of reads for a barcode to be considered (150)
#' @return a ggplot2 histogram or NULL
#'
#' @examples
#' seqsumFile <- system.file('extdata', 'sequencing_summary.txt.bz2', package = 'nanopoRe')
#' importSequencingSummary(seqsumFile)
#' barcodeFile <- system.file('extdata', 'barcoding_summary.txt.bz2', package = 'nanopoRe')
#' SequencingSummaryBarcodeMerge(barcodeFile=barcodeFile)
#' SequenceSummaryBarcodeHistogram()
#'
#' @export
SequenceSummaryBarcodeHistogram <- function(seqsum=NA, bcthreshold = 150) {
    seqsum <- handleBarcodedCache(seqsum)

    barcodes <- 0
    barcodeHistogram <- NULL

    if ("barcode_arrangement" %in% names(seqsum)) {
        barcodedata = plyr::count(seqsum$barcode_arrangement)
        barcodedata = subset(barcodedata, barcodedata$freq > bcthreshold)
        names(barcodedata) <- gsub("x", "barcode", names(barcodedata))
        if ("unclassified" %in% barcodedata$barcode) {
            barcodes <- nrow(barcodedata[-which(barcodedata$barcode == "unclassified"), ])
            barcodeUnass <- sum(barcodedata[-which(barcodedata$barcode == "unclassified"), "freq"])/sum(barcodedata$freq) *
                100
            barcodeRange <- range(subset(barcodedata, "barcode" != "unclassified")$freq)
        } else {
            barcodes <- nrow(barcodedata)
            barcodeUnass = 100
            barcodeRange <- range(barcodedata$freq)
        }
    }

    if (barcodes > 0) {

        barcodeHistogram <- ggplot(barcodedata, aes_string("barcode", "freq", fill = "barcode")) + geom_bar(stat = "identity",
            width = 0.5, fill = "#9ecae1") + xlab("\nDemultiplexed barcodes") + ylab("\nFrequency") +
            scale_y_continuous(expand = c(0, 0)) + labs(title = "Histogram showing abundance of different barcodes") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

    return(ggplot2handler(barcodeHistogram))
}


handleBarcodedCache <- function(seqsum) {
    if (!is.data.frame(seqsum) && is.na(seqsum)) {
        oname <- "barcodedata"
        if (hasCachedObject(oname)) {
            seqsum <- getCachedObject(oname)
        }
    }
    return(seqsum)
}


