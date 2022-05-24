test_that("Tx2Gene", {
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
    object <- Tx2Gene(object)
    expect_s4_class(object, "Tx2Gene")
    con <- file.path(tempdir(), "tx2gene.csv")
    unlink(con, recursive = FALSE)
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
    unlink(con, recursive = FALSE)
})
