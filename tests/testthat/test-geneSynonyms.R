context("geneSynonyms")

test_that("geneSynonyms", {
    organisms <- c("Homo sapiens", "Mus musculus")
    for (organism in organisms) {
        object <- geneSynonyms(organism = organism)
        expect_s4_class(object, "DataFrame")
    }
})
