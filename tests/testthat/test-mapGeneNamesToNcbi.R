test_that("Homo sapiens", {
    expect_identical(
        object = mapGeneNamesToNcbi(
            genes = c("TUT4", "ZCCHC11", "TENT3A"),
            organism = "Homo sapiens"
        ),
        expected = rep(23318L, times = 3L)
    )
})

test_that("Mus musculus", {
    expect_identical(
        object = mapGeneNamesToNcbi(
            genes = c("Nfe2l2", "Nrf2"),
            organism = "Mus musculus"
        ),
        expected = rep(18024L, times = 2L)
    )
})
