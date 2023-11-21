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
    df <- EnsemblToNcbi(
        object = object,
        organism = "Homo sapiens",
        strict = FALSE
    )
    ## FIXME This is failing, need to rethink.
    expect_identical(object, df[["ensemblGeneId"]])
    expect_false(anyNA(df[["ensemblGeneId"]]))
    ## FIXME This is failing, need to rethink.
    expect_true(anyNA(df[["ncbiGeneId"]]))
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
        strict = FALSE
    )
    expect_identical(object, df[["ncbiGeneId"]])
    expect_false(anyNA(df[["ncbiGeneId"]]))
    expect_true(anyNA(df[["ensemblGeneId"]]))

})
