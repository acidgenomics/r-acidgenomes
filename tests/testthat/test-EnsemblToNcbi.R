formats <- eval(formals(`EnsemblToNcbi,character`)[["format"]])

test_that("EnsemblToNcbi : character", {
    genes <- c(
        "ENSG00000063587",
        "ENSG00000004866",
        "ENSG00000000005",
        "ENSG00000000003"
    )
    ## 1:1 mapping.
    object <- EnsemblToNcbi(genes, format = "1:1")
    expect_identical(metadata(object)[["format"]], "1:1")
    expected <- DataFrame(
        "ensemblGeneId" = genes,
        "ncbiGeneId" = c(10838L, 7982L, 64102L, 7105L)
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
    ## Long format (non-unique).
    object <- EnsemblToNcbi(genes, format = "long")
    expect_identical(metadata(object)[["format"]], "long")
    expected <- DataFrame(
        "ensemblGeneId" = c(
            "ENSG00000000003",
            "ENSG00000000005",
            "ENSG00000004866",
            "ENSG00000004866",
            "ENSG00000063587",
            "ENSG00000063587"
        ),
        "ncbiGeneId" = c(
            7105L,
            64102L,
            7982L,
            93655L,
            10838L,
            105373378L
        )
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
})

test_that("NcbiToEnsembl : character", {
    ## These are from the Ensembl return above.
    genes <- c(
        105373378L,
        93655L,
        64102L,
        10838L,
        7982L,
        7105L
    )
    ## 1:1 mapping (of input keys, note the expected Ensembl dupes here).
    ## FIXME Improve coverage here to ensure we're matching ambiguous dupes
    ## correctly from HGNC.
    object <- NcbiToEnsembl(
        object = genes,
        organism = "Homo sapiens",
        format = "1:1"
    )
    expect_identical(metadata(object)[["format"]], "1:1")
    expected <- DataFrame(
        "ncbiGeneId" = genes,
        "ensemblGeneId" = c(
            "ENSG00000063587",
            "ENSG00000004866",
            "ENSG00000000005",
            "ENSG00000063587",
            "ENSG00000004866",
            "ENSG00000000003"
        )
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
    ## Long format (non-unique).
    object <- NcbiToEnsembl(
        object = genes,
        organism = "Homo sapiens",
        format = "long"
    )
    expect_identical(metadata(object)[["format"]], "long")
    expected <- DataFrame(
        "ncbiGeneId" = c(
            7105L,
            7982L,
            10838L,
            64102L,
            93655L,
            105373378L
        ),
        "ensemblGeneId" = c(
            "ENSG00000000003",
            "ENSG00000004866",
            "ENSG00000063587",
            "ENSG00000000005",
            "ENSG00000004866",
            "ENSG00000063587"
        )
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
})
