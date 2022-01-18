context("mapGencodeToEnsembl")

test_that("mapGencodeToEnsembl", {
    expect_identical(
        object = mapGencodeToEnsembl(39L),
        expected = 105L
    )
})
