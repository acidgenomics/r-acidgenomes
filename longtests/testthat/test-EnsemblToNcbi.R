## See also:
## - https://bioinformatics.stackexchange.com/questions/21/

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

## FIXME 151 mismatches between HGNC and Ensembl currently:
## - ENSG00000111215
## - ENSG00000169627
## - ENSG00000170667
## ...

test_that("Mapping consistency between Ensembl and HGNC", {
    x <- EnsemblToNcbi(hgnc)
    x <- as(x, "DFrame")
    metadata(x) <- list()
    rownames(x) <- NULL
    y <- EnsemblToNcbi(ens)
    y <- as(y, "DFrame")
    metadata(y) <- list()
    rownames(y) <- NULL
    genes <- sort(intersect(x[[1L]], y[[1L]]))
    x <- x[match(genes, table = x[[1L]]), ]
    y <- y[match(genes, table = y[[1L]]), ]
    expect_identical(x, y)
    ## > y[which(x[[2]] != y[[2]]), ]
})
