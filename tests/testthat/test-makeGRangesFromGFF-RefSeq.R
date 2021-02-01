context("makeGRangesFromGFF : RefSeq")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint
file <- file.path("cache", "refseq.gff3")

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "RefSeqGenes")
    expect_identical(length(object), 62L)
    expect_identical(names(object)[[1L]], "CICP27")
    expect_identical(
        object = lapply(mcols(object[[1L]]), class),
        expected = list(
            "broadClass" = Rle,
            "description" = Rle,
            "gbkey" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "geneSynonym" = "CompressedCharacterList",
            "pseudo" = Rle,
            "source" = Rle,
            "type" = Rle
        )
    )
})

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "RefSeqTranscripts")
    expect_identical(length(object), 100L)
    expect_identical(names(object)[[1L]], "NM_001005221.2")
    AsIs <- "list"  # nolint
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "exception" = Rle,
            "gbkey" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "geneSynonym" = AsIs,
            "inference" = Rle,
            "modelEvidence" = Rle,
            "product" = Rle,
            "source" = Rle,
            "txId" = Rle,
            "type" = Rle
        )
    )
})
