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
    genes <- c(64102L, 10838L, 7982L, 7105L)
    object <- NcbiToEnsembl(
        object = genes,
        organism = "Homo sapiens",
        format = "1:1"
    )
    expect_identical(metadata(object)[["format"]], "1:1")
    expected <- DataFrame(
        "ncbiGeneId" = genes,
        "ensemblGeneId" = c(
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
            64102L
        ),
        "ensemblGeneId" = c(
            "ENSG00000000003",
            "ENSG00000004866",
            "ENSG00000063587",
            "ENSG00000000005"
        )
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
})

test_that("Invalid keys", {
    expect_error(
        object = EnsemblToNcbi(
            object = "ENSG00000000000",
            organism = "Homo sapiens"
        ),
        regexp = "ENSEMBL"
    )
    expect_error(
        object = NcbiToEnsembl(
            object = 0L,
            organism = "Homo sapiens"
        ),
        regexp = "ENTREZID"
    )
})
