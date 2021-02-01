context("makeGRangesFromGFF : UCSC")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint
file <- gffs[["ucsc_hg38_ensgene_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "UCSCGenes")
    expect_identical(length(object), 64252L)
    expect_identical(names(object)[[1L]], "ENSG00000223972")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "geneId" = Rle,
            "geneName" = Rle
        )
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "UCSCTranscripts")
    expect_identical(length(object), 208239L)
    expect_identical(names(object)[[1L]], "ENST00000456328")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "geneId" = Rle,
            "geneName" = Rle,
            "txBiotype" = Rle,
            "txChrom" = Rle,
            "txEnd" = Rle,
            "txId" = Rle,
            "txName" = Rle,
            "txNumber" = Rle,
            "txStart" = Rle,
            "txStrand" = Rle
        )
    )
})
