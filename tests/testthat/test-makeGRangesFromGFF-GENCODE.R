context("makeGRangesFromGFF : GENCODE")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint
file <- file.path("cache", "gencode.gtf")

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "GencodeGenes")
    expect_identical(length(object), 60L)
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
    expect_identical(length(object), 167L)
    expect_identical(names(object)[[1L]], "ENST00000456328.2")
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

file <- file.path("cache", "gencode.gff3")

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeGenes")
    expect_identical(length(object), 60L)
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
    expect_identical(length(object), 167L)
    expect_identical(names(object)[[1L]], "ENST00000456328.2")
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
