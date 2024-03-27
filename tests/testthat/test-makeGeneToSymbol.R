## FIXME These functions have been removed, consider reworking.

test_that("EnsDb : unmodified", {
    skip_if_not_installed("EnsDb.Hsapiens.v75")
    object <- makeGeneToSymbolFromEnsDb(
        object = "EnsDb.Hsapiens.v75",
        format = "unmodified",
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "GeneToSymbol")
    expect_identical(nrow(object), 63677L)
    expect_identical(
        object = as.data.frame(object)[
            c(
                "ENSG00000225989",
                "ENSG00000238602",
                "ENSG00000254144",
                "ENSG00000273195"
            ),
        ],
        expected = data.frame(
            "geneId" = c(
                "ENSG00000225989",
                "ENSG00000238602",
                "ENSG00000254144",
                "ENSG00000273195"
            ),
            "geneName" = c(
                "ABCF1",
                "5S_rRNA",
                "7SK",
                "ABCA11P"
            ),
            row.names = c(
                "ENSG00000225989",
                "ENSG00000238602",
                "ENSG00000254144",
                "ENSG00000273195"
            )
        )
    )
})

test_that("EnsDb : makeUnique", {
    skip_if_not_installed("EnsDb.Hsapiens.v75")
    object <- makeGeneToSymbolFromEnsDb(
        object = "EnsDb.Hsapiens.v75",
        format = "makeUnique",
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "GeneToSymbol")
    expect_identical(nrow(object), 63677L)
    expect_identical(
        object = as.data.frame(object)[
            c(
                "ENSG00000225989",
                "ENSG00000238602",
                "ENSG00000254144",
                "ENSG00000273195"
            ),
        ],
        expected = data.frame(
            "geneId" = c(
                "ENSG00000225989",
                "ENSG00000238602",
                "ENSG00000254144",
                "ENSG00000273195"
            ),
            "geneName" = c(
                "ABCF1.2",
                "5S_rRNA.2",
                "7SK.2",
                "ABCA11P.1"
            ),
            row.names = c(
                "ENSG00000225989",
                "ENSG00000238602",
                "ENSG00000254144",
                "ENSG00000273195"
            )
        )
    )
})

test_that("EnsDb : 1:1", {
    skip_if_not_installed("EnsDb.Hsapiens.v75")
    object <- makeGeneToSymbolFromEnsDb(
        object = "EnsDb.Hsapiens.v75",
        format = "1:1",
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "GeneToSymbol")
    expect_identical(nrow(object), 56638L)
    expect_identical(
        object = as.data.frame(object)[
            c(
                "ENSG00000197953",
                "ENSG00000251595",
                "ENSG00000260053",
                "ENSG00000073734"
            ),
        ],
        expected = data.frame(
            "geneId" = c(
                "ENSG00000197953",
                "ENSG00000251595",
                "ENSG00000260053",
                "ENSG00000073734"
            ),
            "geneName" = c(
                "AADACL2",
                "ABCA11P",
                "ABCB10P4",
                "ABCB11"
            ),
            row.names = c(
                "ENSG00000197953",
                "ENSG00000251595",
                "ENSG00000260053",
                "ENSG00000073734"
            )
        )
    )
})
