test_that("Homo sapiens : long", {
    object <- goTermsPerGeneName(
        organism = "Homo sapiens",
        geneNames = c("NFE2L2", "MTOR"),
        format = "long"
    )
    expect_s4_class(object, "DFrame")
    expect_null(rownames(object))
    ## FIXME Check that the input returns sorted as expected.
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
    object <- goTermsPerGeneName(
        organism = "Homo sapiens",
        geneNames = c("NFE2L2", "MTOR"),
        format = "split"
    )
    expect_s4_class(object, "SplitDFrameList")
    expect_length(object, 2L)
    expect_named(object, c("NFE2L2", "MTOR"))
    expect_identical(
        object = as.data.frame(object[["MTOR"]][1L, ]),
        expected = data.frame(
            "geneName" = "MTOR",
            "goCategory" = "BP",
            "goId" = "GO:0001558",
            "goName" = "regulation of cell growth"
        )
    )
})

test_that("Homo sapiens : nested", {
    object <- goTermsPerGeneName(
        organism = "Homo sapiens",
        geneNames = c("NFE2L2", "MTOR"),
        format = "nested"
    )
    expect_identical(nrow(object), 2L)
    expect_null(rownames(object))
})
