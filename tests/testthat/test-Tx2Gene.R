## FIXME Need to improve code coverage of NA value handling.



context("Tx2Gene")

test_that("DataFrame", {
    df <- DataFrame(
        "txId" = c(
            "ENST00000635602.1",
            "ENST00000635506.1"
        ),
        "geneId" = c(
            "ENSG00000283061.1",
            "ENSG00000283061.1"
        )
    )
    object <- Tx2Gene(df)
    expect_s4_class(object, "Tx2Gene")
    expect_output(
        object = summary(object),
        regexp = "genes: 1"
    )
})
