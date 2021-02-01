context("makeGRangesFromGFF : UCSC")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint
file <- file.path("cache", "ucsc.gtf")

test_that("GTF genes", {
    stop("FIXME")
})

test_that("GTF transcripts", {
    stop("FIXME")
})
