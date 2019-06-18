context("Invalid/Corrupt Fastq")

test_that("frankenFastq.fastq parses correctly", {
  expect_equal(fastqValidator(system.file("extdata", "frankenFastq.fastq", package = "nanopoRe", mustWork = TRUE)), FALSE)
  expect_equal(getFastqCount(), 608)
  expect_equal(getFastqBases(), 655281)
  expect_equal(getFastqPlusErrorCount(), 6)
  expect_equal(getMalformedFastqHeaderCount(), 7)
  expect_equal(getZeroLengthSequenceCount(), 0)
  expect_equal(getSequenceQualityMismatchCount(), 1)
  expect_equal(getSkippedLineCount(), 19)
})


test_that("frankenFastq.fastq.gz can be fixed", {
  expect_equal(fastqValidator(nanopoRe::fixFastq(system.file("extdata", "frankenFastq.fastq.gz", package = "nanopoRe", mustWork = TRUE), tempfile(pattern="fixfastq_", fileext=".fq"))), TRUE)

  tfile <- tempfile(pattern="fixfastq_", fileext=".fq.gz", tmpdir = tempdir(check=TRUE))
  fixFastq(system.file("extdata", "frankenFastq.fastq.gz", package = "nanopoRe", mustWork = TRUE), tfile)

  expect_equal(fastqValidator(tfile), TRUE)
})
