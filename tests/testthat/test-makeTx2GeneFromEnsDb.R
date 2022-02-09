context("makeTx2GeneFromEnsDb")

skip_if_not_installed("EnsDb.Hsapiens.v75")

test_that("makeTx2GeneFromEnsDb", {
    object <- makeTx2GeneFromEnsDb("EnsDb.Hsapiens.v75")
    expect_s4_class(object, "Tx2Gene")
})
