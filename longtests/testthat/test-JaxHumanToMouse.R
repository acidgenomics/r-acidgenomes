test_that("1:1 unique mapping", {
    object <- JaxHumanToMouse(unique = TRUE)
    expect_s4_class(object, "JaxHumanToMouse")
    expect_identical(nrow(object), 16983L)
    expect_true(hasNoDuplicates(object[["humanGeneName"]]))
    expect_true(hasNoDuplicates(object[["mouseGeneName"]]))
})

test_that("Allow duplicates", {
    object <- JaxHumanToMouse(unique = FALSE)
    expect_s4_class(object, "JaxHumanToMouse")
    expect_identical(nrow(object), 24601L)
    expect_false(hasNoDuplicates(object[["humanGeneName"]]))
    expect_false(hasNoDuplicates(object[["mouseGeneName"]]))
})
