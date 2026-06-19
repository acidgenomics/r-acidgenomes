test_that("mapGoTerms", {
    skip_if_not_installed("ontologyIndex")
    object <- mapGoTerms()
    expect_s4_class(object, "DFrame")
    expect_identical(
        object = as.data.frame(stringsAsFactors = FALSE, object[1L, ]),
        expected = data.frame(
            stringsAsFactors = FALSE,
            id = "GO:0000001",
            name = "mitochondrion inheritance"
        )
    )
})
