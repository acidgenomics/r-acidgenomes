context("export")

test_that("Tx2Gene", {
    object <- DataFrame(
        transcriptID = c(
            "tx0001",
            "tx0002",
            "tx0003",
            "tx0004"
        ),
        geneID = c(
            "gene0001",
            "gene0001",
            "gene0002",
            "gene0002"
        )
    )
    file <- "tx2gene.csv"
    unlink(file)
    expect_false(file.exists(file))
    export(object, file = file)
    expect_true(file.exists(file))
    unlink(file)
})
