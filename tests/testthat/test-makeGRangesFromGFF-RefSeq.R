context("makeGRangesFromGFF : RefSeq")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint
file <- gffs[["refseq_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "RefSeqGenes")
    expect_identical(length(object), 54651L)
    expect_identical(names(object)[[1L]], "A1BG")
    expect_identical(
        object = lapply(mcols(object[[1L]]), class),
        expected = list(
            "broadClass" = Rle,
            "description" = Rle,
            "endRange" = "CompressedCharacterList",
            "exception" = Rle,
            "gbkey" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "geneSynonym" = "CompressedCharacterList",
            "partial" = Rle,
            "pseudo" = Rle,
            "source" = Rle,
            "startRange" = "CompressedCharacterList",
            "translExcept" = "CompressedCharacterList",
            "type" = Rle
        )
    )
})

## FIXME THIS ISNT RETURNING BROAD CLASS, WHICH WE WANT...

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "RefSeqTranscripts")
    expect_identical(length(object), 163975L)
    expect_identical(names(object)[[1L]], "NM_000014.6")
    AsIs <- "list"  # nolint
    expect_identical(
        object = lapply(mcols(object[[1L]]), class),
        expected = list(
            "broadClass" = Rle,
            "description" = Rle,
            "endRange" = AsIs,
            "exception" = Rle,
            "gbkey" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "geneSynonym" = AsIs,
            "inference" = Rle,
            "modelEvidence" = Rle,
            "partial" = Rle,
            "product" = Rle,
            "pseudo" = Rle,
            "source" = Rle,
            "startRange" = AsIs,
            "tag" = Rle,
            "translExcept" = AsIs,
            "txId" = Rle,
            "txName" = Rle,
            "type" = Rle
        )
    )
})
