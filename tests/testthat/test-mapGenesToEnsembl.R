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

## FIXME This test is now failing due to missing URL:
## https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz
## Seems like the NCBI server is currently down.

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
