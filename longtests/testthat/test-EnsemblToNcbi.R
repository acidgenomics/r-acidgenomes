hgnc <- Hgnc()
ens <- makeGRangesFromEnsembl(organism = "Homo sapiens")

test_that("EnsemblToNcbi : all genes", {
    object <- sort(
        x = unique(na.omit(hgnc[["ensemblGeneId"]])),
        decreasing = TRUE
    )
    expect_error(
        object = EnsemblToNcbi(
            object = object,
            organism = "Homo sapiens",
            strict = TRUE
        ),
        regexp = "match failures"
    )
    df <- EnsemblToNcbi(
        object = object,
        organism = "Homo sapiens",
        format = "1:1",
        strict = FALSE
    )
    expect_true(hasNoDuplicates(df[["ensemblGeneId"]]))
    expect_true(identical(df[["ensemblGeneId"]], object))
    expect_false(anyNA(df[["ensemblGeneId"]]))
    expect_true(anyNA(df[["ncbiGeneId"]]))
    df <- EnsemblToNcbi(
        object = object,
        organism = "Homo sapiens",
        format = "long",
        strict = FALSE
    )
    expect_true(hasDuplicates(df[["ensemblGeneId"]]))
    expect_true(areSetEqual(object, df[["ensemblGeneId"]]))
    expect_false(anyNA(df[["ensemblGeneId"]]))
    expect_true(anyNA(df[["ncbiGeneId"]]))
})

test_that("NcbiToEnsembl : all genes", {
    object <- sort(
        x = unique(na.omit(hgnc[["ncbiGeneId"]])),
        decreasing = TRUE
    )
    expect_error(
        object = NcbiToEnsembl(
            object = object,
            organism = "Homo sapiens",
            strict = TRUE
        ),
        regexp = "match failures"
    )
    df <- NcbiToEnsembl(
        object = object,
        organism = "Homo sapiens",
        format = "1:1",
        strict = FALSE
    )
    expect_true(hasNoDuplicates(df[["ncbiGeneId"]]))
    expect_true(identical(df[["ncbiGeneId"]], object))
    expect_false(anyNA(df[["ncbiGeneId"]]))
    expect_true(anyNA(df[["ensemblGeneId"]]))
    df <- NcbiToEnsembl(
        object = object,
        organism = "Homo sapiens",
        format = "long",
        strict = FALSE
    )
    expect_true(hasDuplicates(df[["ncbiGeneId"]]))
    expect_true(areSetEqual(object, df[["ncbiGeneId"]]))
    expect_false(anyNA(df[["ncbiGeneId"]]))
    expect_true(anyNA(df[["ensemblGeneId"]]))
})

test_that("1:1 mapping consistent with HGNC", {
    ## FIXME Add support for Hgnc class EnsemblToNcbi.
    ## FIXME Add support for Hgnc class NcbiToEnsembl.
    ##
    ## FIXME Check that if we take makeGRangesFromEnsembl, create 1:1 mapping
    ## for EnsemblToEntrez, that it matches our expectation.
    ##
    ## FIXME We also may be able to suport EnsemblToNcbi from Our NcbiGeneInfo object.
})
