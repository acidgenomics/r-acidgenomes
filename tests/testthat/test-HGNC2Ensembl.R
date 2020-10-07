context("HGNC2Ensembl")

skip_if_not(hasInternet())

test_that("HGNC2Ensembl", {
    object <- HGNC2Ensembl()
    expect_s4_class(object, "HGNC2Ensembl")
})
