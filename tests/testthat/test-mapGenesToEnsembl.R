test_that("Homo sapiens", {
    expect_identical(
        object = mapGeneNamesToEnsembl(
            genes = c("TUT4", "ZCCHC11", "TENT3A"),
            organism = "Homo sapiens",
            genomeBuild = "GRCh38",
            release = 109L
        ),
        expected = rep("ENSG00000134744", times = 3L)
    )
})

test_that("Mus musculus", {
    expect_identical(
        object = mapGeneNamesToEnsembl(
            genes = c("Nfe2l2", "Nrf2"),
            organism = "Mus musculus",
            genomeBuild = "GRCm39",
            release = 109L
        ),
        expected = rep("ENSMUSG00000015839", times = 2L)
    )
})
