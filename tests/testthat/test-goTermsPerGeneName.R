test_that("Homo sapiens : long", {
    skip_if_not_installed("ontologyIndex")
    object <- goTermsPerGeneName(
        organism = "Homo sapiens",
        geneNames = c("NFE2L2", "MTOR"),
        format = "long"
    )
    expect_s4_class(object, "DFrame")
    expect_null(rownames(object))
    expect_identical(
        object = as.data.frame(stringsAsFactors = FALSE, object[1L, ]),
        expected = data.frame(
            stringsAsFactors = FALSE,
            geneName = "NFE2L2",
            goCategory = "BP",
            goId = "GO:0002931",
            goName = "response to ischemia"
        )
    )
})

test_that("Homo sapiens : split", {
    skip_if_not_installed("ontologyIndex")
    object <- goTermsPerGeneName(
        organism = "Homo sapiens",
        geneNames = c("NFE2L2", "MTOR"),
        format = "split"
    )
    expect_s4_class(object, "SplitDFrameList")
    expect_length(object, 2L)
    expect_named(object, c("NFE2L2", "MTOR"))
    expect_identical(
        object = as.data.frame(stringsAsFactors = FALSE, object[[1L]][1L, ]),
        expected = data.frame(
            stringsAsFactors = FALSE,
            geneName = "NFE2L2",
            goCategory = "BP",
            goId = "GO:0002931",
            goName = "response to ischemia"
        )
    )
})

test_that("Homo sapiens : nested", {
    skip_if_not_installed("ontologyIndex")
    object <- goTermsPerGeneName(
        organism = "Homo sapiens",
        geneNames = c("NFE2L2", "MTOR"),
        format = "nested"
    )
    expect_identical(nrow(object), 2L)
    expect_null(rownames(object))
    expect_identical(
        object = object[["geneName"]],
        expected = c("NFE2L2", "MTOR")
    )
    expect_identical(
        object = object[1L, "goBp"][[1L]][[1L]],
        expected = "GO:0002931 response to ischemia"
    )
    expect_identical(
        object = object[1L, "goCc"][[1L]][[1L]],
        expected = "GO:0000785 chromatin"
    )
    expect_identical(
        object = object[1L, "goMf"][[1L]][[1L]],
        expected = "GO:0000976 transcription cis-regulatory region binding"
    )
})
