context("md5 checksum validation")

test_that("example.fastq.gz has a file checksum", {
  expect_equal(md5sum(system.file("extdata", "example.fastq.gz", package = "nanopoRe")), "78657ac508334d758b0bf959eb366d99")
  expect_equal(md5sum(system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe")), "07272ed62b1425e8b5ce868c05d6e6f4")
})
