context("extra | geneSynonyms")

organisms <- eval(formals("geneSynonyms")[["organism"]])
test_that("geneSynonyms", {
    for (object in organisms) {
        object <- geneSynonyms(organism = organism)
        expect_s4_class(object, "DataFrame")
    }
})
