test_that("Homo sapiens", {
    expect_identical(
        object = mapGeneNamesToHgnc(genes = c("TUT4", "ZCCHC11", "TENT3A")),
        expected = rep(28981L, times = 3L)
    )
})
