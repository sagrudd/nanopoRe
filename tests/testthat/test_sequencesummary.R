context("Oxford Nanopore sequence_summary.txt parsing")

init()

seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
seqsum <- importSequencingSummary(seqsumFile)

test_that("sequencing_summary.txt.bz2 can be parsed", {
  expect_equal(nrow(seqsum), 10000)

  p <- sequencingSummaryPassGauge()

  expect_is(p,"ggplot")

  platform <- sequencingSummaryGetPlatform(seqsum)
  expect_equal(platform, "MinION")
})
