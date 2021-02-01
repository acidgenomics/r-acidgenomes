context("makeGRangesFromEnsembl")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint

test_that("Genes", {
    object <- makeGRangesFromEnsembl(
        organism = organism,
        level = "genes",
        release = ensemblRelease,
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "EnsemblGenes")
    expect_identical(
        object = metadata(object)[c(
            "annotationHubId",
            "genomeBuild",
            "organism",
            "release"
        )],
        expected = list(
            "annotationHubId" = "AH89180",
            "genomeBuild" = "GRCh38",
            "organism" = "Homo sapiens",
            "release" = 102L
        )
    )
    expect_identical(length(object), 68001L)
    expect_identical(
        object = head(names(object), 3L),
        expected = c(
            "ENSG00000228572",
            "ENSG00000182378",
            "ENSG00000226179"
        )
    )
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "canonicalTranscript" = Rle,
            "description" = Rle,
            "entrezId" = "list",
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneIdVersion" = Rle,
            "geneName" = Rle,
            "seqCoordSystem" = Rle
        )
    )
})

## Transcript verion metadata isn't saved in older EnsDb releases, such as v87.
test_that("Transcripts", {
    object <- makeGRangesFromEnsembl(
        organism = organism,
        level = "transcripts",
        release = ensemblRelease,
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblTranscripts")
    expect_identical(length(object), 254853L)
    expect_identical(
        object = head(names(object), n = 2L),
        expected = c("ENST00000635602.1", "ENST00000635506.1")
    )
    expect_identical(
        object = as.character(head(mcols(object)[["geneId"]], n = 2L)),
        expected = c("ENSG00000283061.1", "ENSG00000283061.1")
    )
    expect_identical(
        object = as.character(head(mcols(object)[["geneIdNoVersion"]], n = 2L)),
        expected = c("ENSG00000283061", "ENSG00000283061")
    )
    expect_identical(
        object = as.character(head(mcols(object)[["geneIdVersion"]], n = 2L)),
        expected = c("ENSG00000283061.1", "ENSG00000283061.1")
    )
    expect_identical(
        object = as.character(head(mcols(object)[["txId"]], n = 2L)),
        expected = c("ENST00000635602.1", "ENST00000635506.1")
    )
    expect_identical(
        object = as.character(head(mcols(object)[["txIdNoVersion"]], n = 2L)),
        expected = c("ENST00000635602", "ENST00000635506")
    )
    expect_identical(
        object = as.character(head(mcols(object)[["txIdVersion"]], n = 2L)),
        expected = c("ENST00000635602.1", "ENST00000635506.1")
    )
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "canonicalTranscript" = Rle,
            "description" = Rle,
            "entrezId" = "list",
            "gcContent" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneIdNoVersion" = Rle,
            "geneIdVersion" = Rle,
            "geneName" = Rle,
            "geneSeqEnd" = Rle,
            "geneSeqStart" = Rle,
            "seqCoordSystem" = Rle,
            "txBiotype" = Rle,
            "txCdsSeqEnd" = Rle,
            "txCdsSeqStart" = Rle,
            "txId" = Rle,
            "txIdNoVersion" = Rle,
            "txIdVersion" = Rle,
            "txName" = Rle,
            "txSupportLevel" = Rle
        )
    )
})

test_that("GRCh37 (hg19)", {
    ## Conditionally test if optional EnsDb.Hsapiens.v75 package is installed.
    skip_if_not("EnsDb.Hsapiens.v75" %in% rownames(installed.packages()))
    ## Genes
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "genes",
        genomeBuild = "GRCh37"
    )
    expect_is(object, "GRanges")
    expect_identical(length(object), 64102L)
    expect_identical(head(names(object), 1L), "ENSG00000000003")
    ## Transcripts
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "transcripts",
        genomeBuild = "GRCh37"
    )
    expect_is(object, "GRanges")
    expect_identical(length(object), 215647L)
    expect_identical(head(names(object), 1L), "ENST00000000233")
})

test_that("UCSC identifier matching (hg38)", {
    x <- makeGRangesFromEnsembl(organism = "Homo sapiens", genomeBuild = "hg38")
    expect_s4_class(x, "GRanges")
})

test_that("Organism with 3 words", {
    x <- makeGRangesFromEnsembl(organism = "Canis lupus familiaris")
    expect_s4_class(x, "GRanges")
})

test_that("Invalid parameters", {
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = "Homo sapiens",
            release = 86L
        ),
        regexp = ">= 87"
    )
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = "AAA",
            genomeBuild = "BBB"
        ),
        regexp = "No ID matched on AnnotationHub"
    )
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = c("Homo sapiens", "Mus musculus")
        ),
        regexp = "isString"
    )
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = "Homo sapiens",
            level = "XXX"
        ),
        regexp = "'arg' should be one of \"genes\", \"transcripts\""
    )
})
