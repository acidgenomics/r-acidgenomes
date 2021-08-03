context("makeGRangesFromGFF : bcbio")

test_that("bcbio `ref-transcripts.gtf` file", {
    file <- file.path("cache", "ref-transcripts.gtf")
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "EnsemblGenes")
    expect_length(object, n = 60L)
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "EnsemblTranscripts")
    expect_length(object, n = 167L)
})



context("makeGRangesFromGFF : Ensembl")

skip_if_not(hasInternet())

file <- gffs[["ensembl_grch38_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblGenes")
    expect_identical(length(object), 60675L)
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
            "geneSource" = Rle,
            "source" = Rle,
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
    expect_s4_class(object, "EnsemblTranscripts")
    expect_identical(length(object), 232024L)
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
            "geneSource" = Rle,
            "source" = Rle,
            "tag" = Rle,
            "txBiotype" = Rle,
            "txId" = Rle,
            "txIdNoVersion" = Rle,
            "txIdVersion" = Rle,
            "txName" = Rle,
            "txSource" = Rle,
            "txSupportLevel" = Rle,
            "type" = Rle
        )
    )
})

file <- gffs[["ensembl_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblGenes")
    expect_identical(length(object), 60675L)
    expect_identical(names(object)[[1L]], "ENSG00000223972.5")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "description" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneIdNoVersion" = Rle,
            "geneIdVersion" = Rle,
            "geneName" = Rle,
            "logicName" = Rle,
            "source" = Rle,
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
    expect_s4_class(object, "EnsemblTranscripts")
    expect_identical(length(object), 232024L)
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
            "description" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneIdNoVersion" = Rle,
            "geneIdVersion" = Rle,
            "geneName" = Rle,
            "logicName" = Rle,
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
