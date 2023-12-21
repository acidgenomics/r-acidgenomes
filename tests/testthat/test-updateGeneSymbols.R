test_that("Homo sapiens", {
    object <- updateGeneSymbols(
        geneNames = c("ZCCHC6", "ZCCHC11"),
        organism = "Homo sapiens"
    )
    expect_identical(
        object = object,
        expected = c(
            "ZCCHC6" = "TUT7",
            "ZCCHC11" = "TUT4"
        )
    )
})

test_that("Mus musculus", {
    object <- updateGeneSymbols(
        geneNames = c("Zcchc6", "Zcchc11"),
        organism = "Mus musculus"
    )
    expect_identical(
        object = object,
        expected = c(
            "Zcchc6" = "Tut7",
            "Zcchc11" = "Tut4"
        )
    )
})

test_that("Match failure", {
    expect_error(
        object = updateGeneSymbols(
            geneNames = c("XXX", "YYY"),
            organism = "Homo sapiens"
        ),
        regexp = "Failed to match"
    )
})
