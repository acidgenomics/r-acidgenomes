test_that("geneSynonyms", {
    object <- geneSynonyms(organism = "Homo sapiens")
    expect_s4_class(object, "DataFrame")
    expect_identical(
        object = object["1", "geneSynonyms"],
        expected = CharacterList(c("A1B", "ABG", "GAB", "HYST2477"))
    )
    object <- geneSynonyms(organism = "Mus musculus")
    expect_s4_class(object, "DataFrame")
    expect_identical(
        object = object["11287", "geneSynonyms"],
        expected = CharacterList(c("A1m", "A2m", "MAM"))
    )
})
