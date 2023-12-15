test_that("Homo sapiens : long", {
    object <- goTermsPerGeneName(
        organism = "Homo sapiens",
        format = "long"
    )
    expect_s4_class(object, "DFrame")
    expect_identical(
        object = as.data.frame(
            object[which(object[["geneName"]] == "MTOR")[[1L]], ]
        ),
        expected = data.frame(
            "geneName" = "MTOR",
            "goCategory" = "BP",
            "goId" = "GO:0001558",
            "goName" = "regulation of cell growth"
        )
    )
})

test_that("Homo sapiens : split", {

})

test_that("Homo sapiens : nested", {
})
