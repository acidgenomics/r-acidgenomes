## FIXME Need to improve "geneId" and "geneName" checks here.
## FIXME Check the first element in each return, similar to GFF checks.

context("makeGRangesFromEnsembl")

skip_if_not(hasInternet())

test_that("Genes", {
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "genes",
        release = 102L,
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
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "canonicalTranscript" = "Rle",
            "description" = "Rle",
            "entrezId" = "CompressedIntegerList",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "seqCoordSystem" = "Rle"
        )
    )
})

## Transcript verion metadata isn't saved in older EnsDb releases, such as v87.
test_that("Transcripts", {
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "transcripts",
        release = 102L,
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
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "canonicalTranscript" = "Rle",
            "description" = "Rle",
            "entrezId" = "CompressedIntegerList",
            "gcContent" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdNoVersion" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "geneSeqEnd" = "Rle",
            "geneSeqStart" = "Rle",
            "seqCoordSystem" = "Rle",
            "txBiotype" = "Rle",
            "txCdsSeqEnd" = "Rle",
            "txCdsSeqStart" = "Rle",
            "txId" = "Rle",
            "txIdNoVersion" = "Rle",
            "txIdVersion" = "Rle",
            "txName" = "Rle",
            "txSupportLevel" = "Rle"
        )
    )
})

test_that("GRCh37", {
    ## Conditionally test if optional EnsDb.Hsapiens.v75 package is installed.
    skip_if_not_installed("EnsDb.Hsapiens.v75")
    ## Genes
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "genes",
        genomeBuild = "GRCh37",
        ignoreVersion = TRUE
    )
    expect_is(object, "EnsemblGenes")
    expect_identical(length(object), 64102L)
    expect_identical(head(names(object), 1L), "ENSG00000228572")
    ## Transcripts
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "transcripts",
        genomeBuild = "GRCh37",
        ignoreVersion = TRUE
    )
    expect_is(object, "EnsemblTranscripts")
    expect_identical(length(object), 215647L)
    expect_identical(head(names(object), 1L), "ENST00000478759")
})

test_that("Organism with 3 words", {
    x <- makeGRangesFromEnsembl(
        organism = "Canis lupus familiaris",
        level = "genes"
    )
    expect_s4_class(x, "EnsemblGenes")
})

test_that("Invalid parameters", {
    ## Currently only supports releases back to Ensembl 87.
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = "Homo sapiens",
            release = 86L
        ),
        regexp = "No entry matched on AnnotationHub"
    )
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = "Homo sapiens",
            genomeBuild = "BBB"
        ),
        regexp = "No entry matched on AnnotationHub"
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
