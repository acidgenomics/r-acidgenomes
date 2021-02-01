## FIXME Need to address parsing issues popping up with vroom engine.
## Warning: One or more parsing issues, see `problems()` for details

context("HGNC2Ensembl")

skip_if_not(hasInternet())

test_that("HGNC2Ensembl", {
    object <- HGNC2Ensembl()
    expect_s4_class(object, "HGNC2Ensembl")
})
