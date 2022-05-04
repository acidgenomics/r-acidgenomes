context("makeTx2GeneFromFASTA")

test_that("Ensembl", {
    file <- txFastas[["ensembl"]]
    object <- makeTx2GeneFromFASTA(file)
    expect_s4_class(object, "Tx2Gene")
    expect_identical(nrow(object), 204563L)
})

test_that("FlyBase", {
    file <- txFastas[["flybase"]]
    object <- makeTx2GeneFromFASTA(file)
    expect_s4_class(object, "Tx2Gene")
    expect_identical(nrow(object), 30748L)
})

test_that("GENCODE", {
    file <- txFastas[["gencode"]]
    object <- makeTx2GeneFromFASTA(file)
    expect_s4_class(object, "Tx2Gene")
    expect_identical(nrow(object), 246624L)
})

test_that("WormBase", {
    file <- txFastas[["wormbase"]]
    object <- makeTx2GeneFromFASTA(file)
    expect_s4_class(object, "Tx2Gene")
    expect_identical(nrow(object), 31989L)
})
