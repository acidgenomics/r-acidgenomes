context("mapEnsemblReleaseToURL")

skip_if_not(hasInternet(url = "https://ensembl.org/"))

test_that("Ensembl 100", {
    expect_identical(
        object = mapEnsemblReleaseToURL(100L),
        expected = "https://apr2020.archive.ensembl.org"
    )
})
