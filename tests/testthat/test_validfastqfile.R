context("ValidFastq")

test_that("example.fastq.gz parses correctly", {
  expect_equal(fastqValidator(system.file("extdata", "example.fastq.gz", package = "nanopoRe", mustWork = TRUE)), TRUE)
  expect_equal(getFastqCount(), 625)
  expect_equal(getFastqBases(), 669649)
})
