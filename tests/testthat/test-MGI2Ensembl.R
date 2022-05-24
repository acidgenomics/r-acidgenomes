skip_if_not(hasInternet())

test_that("MGI2Ensembl", {
    object <- MGI2Ensembl()
    expect_s4_class(object, "MGI2Ensembl")
})
