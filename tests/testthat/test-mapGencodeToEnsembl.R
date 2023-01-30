test_that("mapGencodeToEnsembl", {
    expect_identical(
        object = mapGencodeToEnsembl(42L),
        expected = 108L
    )
    expect_identical(
        object = mapGencodeToEnsembl("M31"),
        expected = 108L
    )
})
