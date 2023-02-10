test_that("gencodeReleaseHistory", {
    for (organism in c("Homo sapiens", "Mus musculus")) {
        df <- gencodeReleaseHistory(organism = organism)
        expect_s4_class(df, "DataFrame")
    }
})
