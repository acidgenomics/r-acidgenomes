## This code depends on biomaRt and Ensembl, which can time out.

context("extra | mapEnsemblReleaseToURL")

test_that("NULL input", {
    expect_identical(
        object = mapEnsemblReleaseToURL(NULL),
        expected = "http://useast.ensembl.org"
    )
})

test_that("integer input", {
    expect_identical(
        object = mapEnsemblReleaseToURL(96L),
        expected = "http://apr2019.archive.ensembl.org"
    )
})
