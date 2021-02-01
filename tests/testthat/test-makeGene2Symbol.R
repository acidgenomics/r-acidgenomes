context("makeGene2Symbol")

skip_if_not_installed("EnsDb.Hsapiens.v75")
formats <- eval(formals(makeGene2SymbolFromEnsDb)[["format"]])
test_that("makeGene2SymbolFromEnsDb", {
    for (format in formats) {
        object <- makeGene2SymbolFromEnsDb(
            object = "EnsDb.Hsapiens.v75",
            format = format
        )
        expect_s4_class(object, "Gene2Symbol")
    }
})
rm(formats)
