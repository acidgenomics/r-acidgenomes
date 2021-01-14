context("makeGRangesFromGFF : FlyBase")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint

file <- file.path("cache", "flybase.gtf")

test_that("GTF genes", {
    expect_warning(
        object = makeGRangesFromGFF(file = file, level = "genes"),
        regexp = "Returning without broad class definitions."
    )
    suppressWarnings({
        object <- makeGRangesFromGFF(file = file, level = "genes")
    })
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 39L)
    expect_identical(names(object)[[1L]], "FBgn0000022")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "source" = Rle,
            "type" = Rle
        )
    )
})

test_that("GTF transcripts", {
    expect_warning(
        object = makeGRangesFromGFF(file = file, level = "transcripts"),
        regexp = "Returning without broad class definitions."
    )
    suppressWarnings({
        object <- makeGRangesFromGFF(file = file, level = "transcripts")
    })
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 72L)
    expect_identical(names(object)[[1L]], "FBtr0070000")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "source" = Rle,
            "txId" = Rle,
            "txName" = Rle,
            "type" = Rle
        )
    )
})
