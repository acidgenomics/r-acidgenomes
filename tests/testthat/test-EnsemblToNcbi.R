## FIXME Move everything from longtests here.
## FIXME Need to add testing for ambiguous matches.

test_that("character : 1:1 (unique)", {
    genes <- c(
        "ENSG00000063587.15",
        "ENSG00000004866.22",
        "ENSG00000000005.6",
        "ENSG00000000003.16"
    )
    object <- EnsemblToNcbi(genes)
    expected <- DataFrame(
        "ensemblGeneId" = c(
            "ENSG00000063587",
            "ENSG00000004866",
            "ENSG00000000005",
            "ENSG00000000003"
        ),
        "ncbiGeneId" = c(
            10838L,
            7982L,
            64102L,
            7105L
        ),
        row.names = genes
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
})

test_that("character : Invalid key", {
    expect_error(
        object = EnsemblToNcbi(
            object = "ENSG00000000000",
            organism = "Homo sapiens"
        ),
        regexp = "ENSEMBL"
    )
})
