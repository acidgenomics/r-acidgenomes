context("export")

## FIXME Need to check that this doesn't export with column names.

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
    con <- file.path(tempdir(), "tx2gene.csv")
    unlink(file)
    expect_false(file.exists(con))
    export(object, con = con)
    expect_true(file.exists(con))
    unlink(con)
})
