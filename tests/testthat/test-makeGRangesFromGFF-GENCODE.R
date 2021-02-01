context("makeGRangesFromGFF : GENCODE")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint
file <- gffs[["gencode_grch38_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "GencodeGenes")
    expect_identical(length(object), 60660L)
    expect_identical(names(object)[[1L]], "ENSG00000223972.5")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneIdNoVersion" = Rle,
            "geneIdVersion" = Rle,
            "geneName" = Rle,
            "havanaGene" = Rle,
            "level" = Rle,
            "source" = Rle,
            "tag" = Rle,
            "type" = Rle
        )
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeTranscripts")
    expect_identical(length(object), 232117L)
    expect_identical(
        object = names(object)[[1L]],
        expected = "ENST00000456328.2"
    )
    expect_identical(
        object = as.character(mcols(object)[["txId"]])[[1L]],
        expected = "ENST00000456328.2"
    )
    expect_identical(
        object = as.character(mcols(object)[["txIdNoVersion"]])[[1L]],
        expected = "ENST00000456328"
    )
    expect_identical(
        object = as.character(mcols(object)[["txIdVersion"]])[[1L]],
        expected = "ENST00000456328.2"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneId"]])[[1L]],
        expected = "ENSG00000223972.5"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneIdNoVersion"]])[[1L]],
        expected = "ENSG00000223972"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneIdVersion"]])[[1L]],
        expected = "ENSG00000223972.5"
    )
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "ccdsId" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneIdNoVersion" = Rle,
            "geneIdVersion" = Rle,
            "geneName" = Rle,
            "havanaGene" = Rle,
            "havanaTranscript" = Rle,
            "hgncId" = Rle,
            "level" = Rle,
            "ont" = Rle,
            "proteinId" = Rle,
            "source" = Rle,
            "tag" = Rle,
            "txBiotype" = Rle,
            "txId" = Rle,
            "txIdNoVersion" = Rle,
            "txIdVersion" = Rle,
            "txName" = Rle,
            "txSupportLevel" = Rle,
            "type" = Rle
        )
    )
})

file <- gffs[["gencode_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeGenes")
    expect_identical(length(object), 60660L)
    expect_identical(names(object)[[1L]], "ENSG00000223972.5")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneIdNoVersion" = Rle,
            "geneIdVersion" = Rle,
            "geneName" = Rle,
            "havanaGene" = Rle,
            "hgncId" = Rle,
            "level" = Rle,
            "ont" = "CompressedCharacterList",
            "source" = Rle,
            "tag" = "CompressedCharacterList",
            "type" = Rle
        )
    )
})

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeTranscripts")
    expect_identical(length(object), 232117L)
    expect_identical(
        object = names(object)[[1L]],
        expected = "ENST00000456328.2"
    )
    expect_identical(
        object = as.character(mcols(object)[["txId"]])[[1L]],
        expected = "ENST00000456328.2"
    )
    expect_identical(
        object = as.character(mcols(object)[["txIdNoVersion"]])[[1L]],
        expected = "ENST00000456328"
    )
    expect_identical(
        object = as.character(mcols(object)[["txIdVersion"]])[[1L]],
        expected = "ENST00000456328.2"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneId"]])[[1L]],
        expected = "ENSG00000223972.5"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneIdNoVersion"]])[[1L]],
        expected = "ENSG00000223972"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneIdVersion"]])[[1L]],
        expected = "ENSG00000223972.5"
    )
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "ccdsId" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneIdNoVersion" = Rle,
            "geneIdVersion" = Rle,
            "geneName" = Rle,
            "havanaGene" = Rle,
            "havanaTranscript" = Rle,
            "hgncId" = Rle,
            "level" = Rle,
            "ont" = "CompressedCharacterList",
            "proteinId" = Rle,
            "source" = Rle,
            "tag" = "CompressedCharacterList",
            "txBiotype" = Rle,
            "txId" = Rle,
            "txIdNoVersion" = Rle,
            "txIdVersion" = Rle,
            "txName" = Rle,
            "txSupportLevel" = Rle,
            "type" = Rle
        )
    )
})
