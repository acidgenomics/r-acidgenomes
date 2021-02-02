context("Ensembl2Entrez")

formats <- eval(formals(`Ensembl2Entrez,GRanges`)[["format"]])

test_that("character", {
    genes <- c(
        "ENSG00000000003",
        "ENSG00000000005",
        "ENSG00000004866",
        "ENSG00000063587"
    )
    ## 1:1 mapping.
    object <- Ensembl2Entrez(genes, format = "1:1")
    expect_identical(metadata(object)[["format"]], "1:1")
    expected <- DataFrame(
        "ensemblId" = genes,
        "entrezId" = c(7105L, 64102L, 7982L, 10838L),
        row.names = genes
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
    ## Long format (non-unique).
    object <- Ensembl2Entrez(genes, format = "long")
    expect_identical(metadata(object)[["format"]], "long")
    expected <- DataFrame(
        "ensemblId" = c(
            "ENSG00000000003",
            "ENSG00000000005",
            "ENSG00000004866",
            "ENSG00000004866",
            "ENSG00000063587",
            "ENSG00000063587"
        ),
        "entrezId" = c(
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



context("Entrez2Ensembl")

test_that("character", {
    ## These are from the Ensembl return above.
    genes <- c(7105L, 7982L, 10838L, 64102L, 93655L, 105373378L)
    ## 1:1 mapping (of input keys, note the expected Ensembl dupes here).
    object <- Entrez2Ensembl(
        object = genes,
        organism = "Homo sapiens",
        format = "1:1"
    )
    expect_identical(metadata(object)[["format"]], "1:1")
    expected <- DataFrame(
        "entrezId" = genes,
        "ensemblId" = c(
            "ENSG00000000003",
            "ENSG00000004866",
            "ENSG00000063587",
            "ENSG00000000005",
            "ENSG00000004866",
            "ENSG00000063587"
        ),
        row.names = genes
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
    ## Long format (non-unique).
    object <- Entrez2Ensembl(
        object = genes,
        organism = "Homo sapiens",
        format = "long"
    )
    expect_identical(metadata(object)[["format"]], "long")
    expected <- DataFrame(
        "entrezId" = genes,
        "ensemblId" = c(
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
