test_that("mapGoTerms", {
    object <- mapGoTerms()
    expect_s4_class(object, "DFrame")
    expect_identical(
        object = as.data.frame(object[1L, ]),
        expected = data.frame(
            "id" = "GO:0000001",
            "name" = "mitochondrion inheritance"
        )
    )
})
