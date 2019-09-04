context("Reference Genome handling")

init()
referenceFasta <- system.file("extdata",
                              "Escherichia_coli_complete_genome.fasta",
                              package = "nanopoRe")

test_that("setReferenceGenome -- reference genome can be set", {
    expect_silent(setReferenceGenome(referenceFasta))
})

test_that("getReferenceGenome -- reference genome can be retrieved", {
    setReferenceGenome(referenceFasta)
    expect_true(getReferenceGenome() == referenceFasta)
})

test_that("loadReferenceGenome -- parsing stuff in memory", {
    setReferenceGenome(referenceFasta)
    expect_silent(loadReferenceGenome())
})

test_that("getChromosomeIds -- can extract chromosome names", {
    setReferenceGenome(referenceFasta)
    loadReferenceGenome()
    expect_vector(getChromosomeIds())
    expect_equal(getChromosomeIds(), c("Escherichia_coli_plasmid",
                                       "Escherichia_coli_chromosome"))
})


test_that("getSeqLengths -- getting reference sequence lengths", {
    setReferenceGenome(referenceFasta)
    loadReferenceGenome()
    expect_vector(getSeqLengths(getChromosomeIds()))
    expect_equal(length(getSeqLengths(getChromosomeIds())), 2)
    expect_equal(as.numeric(getSeqLengths(getChromosomeIds())), c(110007,4765434))
})


test_that("chromosomeMappingSummary -- quick summary stats", {
    demoBam <- system.file("extdata", "Ecoli_zymo_R10_filt_subs.bam", package = "nanopoRe")
    setReferenceGenome(referenceFasta)
    loadReferenceGenome()
    mapStats <- chromosomeMappingSummary(getChromosomeIds(), demoBam)

    # chrId chrLength N.... GC.... Mapped.Reads Mapped.Bases Mean.Coverage
    # 1 Escherichia_coli_chromosome 4,765,434     0  50.72          802    3,522,093          0.75
    # 2    Escherichia_coli_plasmid   110,007     0  49.15           18       70,012          0.65
    expect_equal(class(mapStats), "data.frame")
    expect_equal(dim(mapStats), c(2, 7))
})


test_that("getStringSetId -- pointers for reference sequences", {
    setReferenceGenome(referenceFasta)
    loadReferenceGenome()
    ssid <- getStringSetId(getChromosomeIds())
    expect_vector(ssid)
    expect_equal(ssid, c(1, 2))
    expect_equal(getStringSetId("Escherichia_coli_chromosome"), 2)
})


test_that("cleanReferenceGenome -- removing sequence pointers", {
    setReferenceGenome(referenceFasta)
    loadReferenceGenome()
    expect_equal(length(getSeqLengths(getChromosomeIds())), 2)
    expect_silent(cleanReferenceGenome(getChromosomeIds()[1]))
    expect_equal(length(getSeqLengths(getChromosomeIds())), 1)
    expect_equal(getChromosomeIds(), "Escherichia_coli_chromosome")
})

test_that("getChromosomeSequence -- getting the DNA", {
    setReferenceGenome(referenceFasta)
    loadReferenceGenome()
    expect_silent(dna <- getChromosomeSequence(getStringSetId("Escherichia_coli_chromosome")))
    expect_type(dna, "S4")
    expect_equal("DNAString" %in% class(dna), TRUE)
    expect_equal(length(dna), 4765434)
})


