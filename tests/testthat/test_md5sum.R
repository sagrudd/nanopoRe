context("md5 checksum validation")

test_that("example.fastq.gz has a file checksum", {
  expect_equal(nanopoRe::md5sum(system.file("extdata", "example.fastq.gz", package = "nanopoRe", mustWork = TRUE)), "78657ac508334d758b0bf959eb366d99")
})
