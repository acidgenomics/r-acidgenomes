## This code depends on biomaRt and Ensembl, which can time out.

context("extra | matchEnsemblReleaseToURL")

test_that("NULL input", {
    expect_identical(
        object = matchEnsemblReleaseToURL(NULL),
        expected = "http://useast.ensembl.org"
    )
})

test_that("integer input", {
    expect_identical(
        object = matchEnsemblReleaseToURL(96L),
        expected = "http://apr2019.archive.ensembl.org"
    )
})
