context("BAM mapping file parsing and review")

init()
demoBam <- system.file("extdata", "Ecoli_zymo_R10_filt_subs.bam", package = "nanopoRe")
referenceFasta <- system.file("extdata",
                              "Escherichia_coli_complete_genome.fasta",
                              package = "nanopoRe")
setReferenceGenome(referenceFasta)
demoBamIdx <- system.file("extdata",
                          "Ecoli_zymo_R10_filt_subs.bam",
                          package = "nanopoRe")

test_that("testBam -- demo BAM file can be parsed", {
    bamData <- testBam(demoBam, yieldSize = 10L)

    expect_type(bamData, "list")
    expect_equal(names(bamData), c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar", "qual", "tag"))
    expect_type(bamData$qual, "S4")
    expect_equal("PhredQuality" %in% class(bamData$qual), TRUE)
    expect_equal(length(bamData$qual), 10)
})

test_that("bamSummarise -- can prepare mapping summary information", {
    suppressMessages(bamSummary <- bamSummarise(demoBam, force=FALSE, blockSize=1000L))

    expect_equal(class(bamSummary), "data.frame")
    expect_equal(dim(bamSummary), c(867, 12))
    expect_equal(colnames(bamSummary), c("qname", "readFlag", "rname", "strand", "start", "qwidth", "end", "mapq", "readq", "coverage", "accuracy", "identity"))
})

test_that("bamSummariseByChr -- check functionality", {
    suppressMessages(bamSummaryChr <- bamSummariseByChr(chrId=getChromosomeIds()[1],
        bamFile=demoBam,
        blockSize=100L,
        index=demoBamIdx))

    expect_equal(class(bamSummaryChr), "data.frame")
    expect_equal(dim(bamSummaryChr), c(19, 12))
    expect_equal(colnames(bamSummaryChr), c("qname", "readFlag", "rname", "strand", "start", "qwidth", "end", "mapq", "readq", "coverage", "accuracy", "identity"))

})


test_that("bamSummaryToCoverage -- extraction of GRanges", {
    bamGR <- bamSummaryToCoverage(demoBam, tilewidth=1000)

    expect_type(bamGR, "S4")
    expect_equal("GRanges" %in% class(bamGR), TRUE)
    expect_equal(length(bamGR), 4877)
})
