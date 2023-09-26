test_that("makeGeneToSymbolFromEnsDb", {
    skip_if_not_installed("EnsDb.Hsapiens.v75")
    formats <- eval(formals(makeGeneToSymbolFromEnsDb)[["format"]])
    for (format in formats) {
        object <- makeGeneToSymbolFromEnsDb(
            object = "EnsDb.Hsapiens.v75",
            format = format
        )
        expect_s4_class(object, "GeneToSymbol")
    }
})
