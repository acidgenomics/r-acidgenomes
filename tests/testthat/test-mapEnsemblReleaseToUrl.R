## This code depends on Ensembl server, which can time out.

test_that("Ensembl 100", {
    expect_identical(
        object = mapEnsemblReleaseToUrl(100L),
        expected = "https://apr2020.archive.ensembl.org"
    )
})

test_that("NULL input", {
    expect_identical(
        object = mapEnsemblReleaseToUrl(NULL),
        expected = "https://useast.ensembl.org"
    )
})
