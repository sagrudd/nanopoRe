context("ValidFastq")

fastqFile <- system.file("extdata", "example.fastq.gz", package = "nanopoRe")

test_that("example.fastq.gz parses correctly", {
  expect_equal(fastqValidator(fastqFile), TRUE)
  expect_equal(getFastqCount(), 625)
  expect_equal(getFastqBases(), 669649)
})


test_that("fastqCheckup -- looking at fastq file behaviour", {
    expect_silent(fqcheck <- fastqCheckup(fastqFile))
    expect_vector(fqcheck)
    expect_equal(names(fqcheck), c("filename","checksum","reads","bases","size","nt","n" ))
    expect_equal(as.numeric(fqcheck[['reads']]), 625)
})
