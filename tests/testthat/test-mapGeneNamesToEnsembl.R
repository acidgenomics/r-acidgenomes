test_that("Homo sapiens", {
    expect_identical(
        object = mapGeneNamesToEnsembl(
            genes = c("TUT4", "ZCCHC11", "TENT3A"),
            organism = "Homo sapiens"
        ),
        expected = rep("ENSG00000134744", times = 3L)
    )
})

test_that("Mus musculus", {
    expect_identical(
        object = mapGeneNamesToEnsembl(
            genes = c("Nfe2l2", "Nrf2"),
            organism = "Mus musculus"
        ),
        expected = rep("ENSMUSG00000015839", times = 2L)
    )
})
