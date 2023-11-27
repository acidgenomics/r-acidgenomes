## FIXME Move everything from longtests here.
## FIXME Need to add testing for ambiguous matches.

test_that("character : 1:1 (unique)", {
    genes <- c(64102L, 10838L, 7982L, 7105L)
    object <- NcbiToEnsembl(genes, organism = "Homo sapiens")
    expected <- DataFrame(
        "ncbiGeneId" = genes,
        "ensemblGeneId" = c(
            "ENSG00000000005",
            "ENSG00000063587",
            "ENSG00000004866",
            "ENSG00000000003"
        ),
        row.names = genes
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
})

test_that("Invalid key", {
    expect_error(
        object = NcbiToEnsembl(
            object = 0L,
            organism = "Homo sapiens"
        ),
        regexp = "ENTREZID"
    )
})
