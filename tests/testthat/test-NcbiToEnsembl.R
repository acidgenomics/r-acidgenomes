## FIXME Move everything from longtests here.
## FIXME Need to add testing for ambiguous matches.

test_that("character : 1:1 (unique)", {
    object <- NcbiToEnsembl(
        object = c(
            64102L,
            10838L,
            7982L,
            7105L
        ),
        organism = "Homo sapiens",
        format = "1:1"
    )
    expect_identical(metadata(object)[["format"]], "1:1")
    expected <- DataFrame(
        "ncbiGeneId" = c(
            64102L,
            10838L,
            7982L,
            7105L
        ),
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
})

## FIXME Need to use a gene that multi-maps here.

test_that("character : long (non-unique)", {
    object <- NcbiToEnsembl(
        object = c(
            64102L,
            10838L,
            7982L,
            7105L
        ),
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

test_that("Invalid key", {
    expect_error(
        object = NcbiToEnsembl(
            object = 0L,
            organism = "Homo sapiens"
        ),
        regexp = "ENTREZID"
    )
})
