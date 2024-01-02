test_that("Homo sapiens : GRCh38 : exons", {
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "exons",
        release = 110L,
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "EnsemblExons")
    expect_length(
        object = object,
        n = n[["hsapiens"]][["ensembl"]][["exons"]]
    )
    gr <- unlist(object, recursive = FALSE, use.names = FALSE)
    expect_length(
        object = unique(mcols(gr)[["txId"]]),
        n = n[["hsapiens"]][["ensembl"]][["transcripts"]]
    )
    expect_length(
        object = unique(mcols(gr)[["geneId"]]),
        n = n[["hsapiens"]][["ensembl"]][["genes"]]
    )
})

test_that("Homo sapiens : GRCh37 : genes", {
    skip_if_not_installed("EnsDb.Hsapiens.v75")
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "genes",
        genomeBuild = "GRCh37",
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "EnsemblGenes")
    expect_length(object, 64102L) # FIXME
    ## FIXME Rework to use expect_named.
    expect_identical(head(names(object), 1L), "ENSG00000228572")
    ## FIXME Ensure all names match expected pattern.
})

test_that("Homo sapiens : GRCh37 : transcripts", {
    skip_if_not_installed("EnsDb.Hsapiens.v75")
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "transcripts",
        genomeBuild = "GRCh37",
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "EnsemblTranscripts")
    expect_length(object, 215647L) # FIXME
    ## FIXME Rework to use expect_named.
    expect_identical(head(names(object), 1L), "ENST00000478759")
    ## FIXME Ensure all names match expected pattern.
})

## FIXME Ensure we cover exons of GRCh37 as well.

test_that("Mus musculus : GRCm39 : exons", {
    object <- makeGRangesFromEnsembl(
        organism = "Mus musculus",
        level = "exons",
        release = 110L,
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "EnsemblExons")
    expect_length(
        object = object,
        n = n[["mmusculus"]][["ensembl"]][["exons"]]
    )
})

test_that("Legacy GRCh38 87 release without gene version", {
    ## FIXME Check that this errors if user sets ignoreVersion = TRUE.
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "genes",
        genomeBuild = "GRCh38",
        release = 87L,
        ignoreVersion = FALSE,
        extraMcols = TRUE
    )
    expect_s4_class(object, "EnsemblGenes")
})

test_that("Organism with 3 words", {
    x <- makeGRangesFromEnsembl(
        organism = "Canis lupus familiaris",
        level = "genes",
        release = 110L,
        ignoreVersion = FALSE
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
