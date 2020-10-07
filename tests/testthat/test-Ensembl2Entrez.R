context("Ensembl2Entrez")

formats <- eval(formals(`Ensembl2Entrez,GRanges`)[["format"]])

test_that("character", {
    genes <- c(
        "ENSG00000000003", "ENSG00000000005",
        "ENSG00000004866", "ENSG00000063587"
    )
    ## 1:1 mapping.
    object <- Ensembl2Entrez(genes, format = "1:1")
    expect_identical(metadata(object)[["format"]], "1:1")
    expected <- DataFrame(
        ensembl = genes,
        entrez = c(7105L, 64102L, 7982L, 10838L),
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
        ensembl = c(
            "ENSG00000000003",
            "ENSG00000000005",
            "ENSG00000004866", "ENSG00000004866",
            "ENSG00000063587", "ENSG00000063587"
        ),
        entrez = c(
            7105L,
            64102L,
            7982L, 93655L,
            10838L, 105373378L
        )
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
})

test_that("GRanges", {
    for (format in formats) {
        object <- Ensembl2Entrez(object = rowRanges(rse), format = format)
        expect_s4_class(object, "Ensembl2Entrez")
    }
})

test_that("RangedSummarizedExperiment", {
    for (format in formats) {
        object <- Ensembl2Entrez(object = rse, format = format)
        expect_s4_class(object, "Ensembl2Entrez")
    }
})



context("Entrez2Ensembl")

test_that("character", {
    ## These are from the Ensembl return above.
    genes <- c(
        7105L, 7982L, 10838L,
        64102L, 93655L, 105373378L
    )
    ## 1:1 mapping (of input keys, note the expected Ensembl dupes here).
    object <- Entrez2Ensembl(genes, format = "1:1")
    expect_identical(metadata(object)[["format"]], "1:1")
    expected <- DataFrame(
        entrez = genes,
        ensembl = c(
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
    object <- Entrez2Ensembl(genes, format = "long")
    expect_identical(metadata(object)[["format"]], "long")
    expected <- DataFrame(
        entrez = genes,
        ensembl = c(
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
