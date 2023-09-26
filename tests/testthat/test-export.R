test_that("TxToGene", {
    object <- DataFrame(
        "txId" = c(
            "tx0001",
            "tx0002",
            "tx0003",
            "tx0004"
        ),
        "geneId" = c(
            "gene0001",
            "gene0001",
            "gene0002",
            "gene0002"
        )
    )
    object <- TxToGene(object)
    expect_s4_class(object, "TxToGene")
    con <- file.path(tempdir2(), "tx2gene.csv")
    expect_false(file.exists(con))
    export(object = object, con = con)
    expect_true(file.exists(con))
    expect_identical(
        object = import(
            con = con,
            format = "lines",
            nMax = 4L
        ),
        expected = c(
            "\"tx0001\",\"gene0001\"",
            "\"tx0002\",\"gene0001\"",
            "\"tx0003\",\"gene0002\"",
            "\"tx0004\",\"gene0002\""
        )
    )
    unlink2(con)
})
