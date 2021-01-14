context("makeGRangesFromGFF : WormBase")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint

file <- file.path("cache", "wormbase.gtf")

test_that("GTF genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 66L)
    expect_identical(names(object)[[1L]], "WBGene00000674")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneSource" = Rle,
            "source" = Rle,
            "type" = Rle
        )
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 83L)
    expect_identical(names(object)[[1L]], "B0545.1a.1")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneSource" = Rle,
            "source" = Rle,
            "transcriptBiotype" = Rle,
            "transcriptId" = Rle,
            "transcriptSource" = Rle,
            "type" = Rle
        )
    )
})
