context("makeGRangesFromGFF : GENCODE")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint

file <- file.path("cache", "gencode.gtf")

test_that("GTF genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 60L)
    expect_identical(names(object)[[1L]], "ENSG00000177757.2")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
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
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 167L)
    expect_identical(names(object)[[1L]], "ENST00000326734.2")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "ccdsId" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "havanaGene" = Rle,
            "havanaTranscript" = Rle,
            "level" = Rle,
            "ont" = Rle,
            "proteinId" = Rle,
            "source" = Rle,
            "tag" = Rle,
            "transcriptBiotype" = Rle,
            "transcriptId" = Rle,
            "transcriptName" = Rle,
            "transcriptSupportLevel" = Rle,
            "type" = Rle
        )
    )
})

file <- file.path("cache", "gencode.gff3")

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 60L)
    expect_identical(names(object)[[1L]], "ENSG00000177757.2")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "havanaGene" = Rle,
            "level" = Rle,
            "source" = Rle,
            "type" = Rle
        )
    )
})

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 167L)
    expect_identical(names(object)[[1L]], "ENST00000326734.2")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "ccdsId" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "havanaGene" = Rle,
            "havanaTranscript" = Rle,
            "level" = Rle,
            "proteinId" = Rle,
            "source" = Rle,
            "transcriptBiotype" = Rle,
            "transcriptId" = Rle,
            "transcriptName" = Rle,
            "transcriptSupportLevel" = Rle,
            "type" = Rle
        )
    )
})
