test_that("Hgnc", {
    object <- Hgnc()
    expect_s4_class(object, "Hgnc")
})
