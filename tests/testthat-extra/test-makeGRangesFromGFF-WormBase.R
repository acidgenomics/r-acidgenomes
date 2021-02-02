context("extra | makeGRangesFromGFF : WormBase")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint
file <- gffs[["wormbase_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 46934L)
    expect_identical(names(object)[[1L]], "WBGene00022276")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "geneSource" = Rle,
            "source" = Rle,
            "type" = Rle
        )
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 59897L)
    expect_identical(names(object)[[1L]], "Y74C9A.2a.3")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "geneSource" = Rle,
            "source" = Rle,
            "txBiotype" = Rle,
            "txId" = Rle,
            "txName" = Rle,
            "txSource" = Rle,
            "type" = Rle
        )
    )
})
