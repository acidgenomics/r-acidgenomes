context("geneSynonyms")

skip_if_not(hasInternet())

test_that("geneSynonyms", {
    expect_s4_class(
        object = geneSynonyms(organism = "Homo sapiens"),
        class = "DataFrame"
    )
})
