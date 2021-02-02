context("makeGRangesFromGFF : FlyBase")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint
file <- gffs[["flybase_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "FlyBaseGenes")
    expect_identical(length(object), 17875L)
    expect_identical(names(object)[[1L]], "FBgn0031208")
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
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "FlyBaseTranscripts")
    expect_identical(length(object), 35643L)
    expect_identical(names(object)[[1L]], "FBtr0475186")
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
