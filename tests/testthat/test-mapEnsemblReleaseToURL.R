## This code depends on Ensembl server, which can time out.

skip_if_not(hasInternet(url = "https://ensembl.org/"))

test_that("Ensembl 100", {
    expect_identical(
        object = mapEnsemblReleaseToURL(100L),
        expected = "https://apr2020.archive.ensembl.org"
    )
})

test_that("NULL input", {
    expect_identical(
        object = mapEnsemblReleaseToURL(NULL),
        expected = "http://useast.ensembl.org"
    )
})
