context("Oxford Nanopore sequence_summary.txt parsing")

init()
library(emojifont)


test_that("sequencing_summary.txt.bz2 can be parsed", {
  seqsumFile <- system.file("extdata", "sequencing_summary.txt.bz2", package = "nanopoRe", mustWork = TRUE)
  seqsum <<- importSequencingSummary(seqsumFile)

  expect_equal(nrow(seqsum), 10000)

  platform <- SequencingSummaryGetPlatform(seqsum)
  expect_equal(platform, "MinION")
})

test_that("passed and failed analytics make sense", {
  passedSeqs <<- seqsum[which(seqsum$passes_filtering), ]
  expect_equal(nrow(passedSeqs), 8355)
  failedSeqs <- seqsum[which(!seqsum$passes_filtering), ]
  expect_equal(nrow(failedSeqs), 1645)
})


test_that("SequencingSummary based plots", {

  plot1 <- SequencingSummaryPassGauge(seqsum)
  expect_is(plot1,"ggplot")

  plot2 <- SequencingSummaryChannelActivity(seqsum)
  expect_is(plot2,"ggplot")

  plot3 <- SequencingSummaryWeightedReadLength(seqsum)
  expect_is(plot3,"ggplot")

  plot4 <- SequencingSummaryReadLengthHistogram(seqsum)
  expect_is(plot4,"ggplot")

  plot5 <- SequencingSummaryReadQualityHistogram(seqsum)
  expect_is(plot5,"ggplot")

  plot6 <- SequencingSummaryReadLengthQualityDensity(seqsum)
  expect_is(plot6,"ggplot")

  plot7 <- SequencingSummaryTemporalThroughput(seqsum)
  expect_is(plot7,"ggplot")

  plot8 <- SequencingSummaryCumulativeBases(seqsum)
  expect_is(plot8,"ggplot")

  plot9 <- SequencingSummaryCumulativeReads(seqsum)
  expect_is(plot9,"ggplot")

  plot10 <- SequencingSummarySpeedPlot(seqsum)
  expect_is(plot10,"ggplot")

  plot11 <- SequencingSummaryActiveChannelPlot(seqsum)
  expect_is(plot11,"ggplot")

  file1 <- SequenceSummaryExecutiveSummary(seqsum)
  expect_equal(file.exists(file1), TRUE)

  file2 <- SequenceSummaryBasicInfoPlot(seqsum)
  expect_equal(file.exists(file2), TRUE)
})
