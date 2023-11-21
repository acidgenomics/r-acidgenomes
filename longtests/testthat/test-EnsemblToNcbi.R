hgnc <- Hgnc()

test_that("EnsemblToNcbi : all genes", {
    object <- sort(unique(na.omit(hgnc[["ensemblGeneId"]])))
    expect_error(
        object = EnsemblToNcbi(
            object = object,
            organism = "Homo sapiens",
            strict = TRUE
        ),
        regexp = "match failures"
    )
    ## FIXME This is failing return...
    df <- EnsemblToNcbi(
        object = object,
        organism = "Homo sapiens",
        format = "1:1",
        strict = FALSE
    )
    ## FIXME This is failing, need to rethink.
    expect_identical(object, df[["ensemblGeneId"]])
    ## FIXME This is now messed argh.
    expect_false(anyNA(df[["ensemblGeneId"]]))
    ## FIXME This is failing, need to rethink.
    expect_true(anyNA(df[["ncbiGeneId"]]))

    ## FIXME Add coverage for long mode.
    df <- EnsemblToNcbi(
        object = object,
        organism = "Homo sapiens",
        format = "long",
        strict = FALSE
    )
})

test_that("NcbiToEnsembl : all genes", {
    object <- sort(unique(na.omit(hgnc[["ncbiGeneId"]])))
    expect_error(
        object = NcbiToEnsembl(
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
    expect_identical(object, df[["ncbiGeneId"]])
    expect_false(anyNA(df[["ncbiGeneId"]]))
    expect_true(anyNA(df[["ensemblGeneId"]]))

    ## FIXME Add coverage for this.
    df <- EnsemblToNcbi(
        object = object,
        organism = "Homo sapiens",
        format = "long",
        strict = FALSE
    )
})
