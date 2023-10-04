test_that("1:1 unique mapping", {
    object <- HumanToMouse(unique = TRUE)
    expect_s4_class(object, "HumanToMouse")
    expect_identical(nrow(object), 16985L)
    expect_true(hasNoDuplicates(object[["humanGeneName"]]))
    expect_true(hasNoDuplicates(object[["mouseGeneName"]]))
})

test_that("Allow duplicates", {
    object <- HumanToMouse(unique = FALSE)
    expect_s4_class(object, "HumanToMouse")
    expect_identical(nrow(object), 24600L)
    expect_false(hasNoDuplicates(object[["humanGeneName"]]))
    expect_false(hasNoDuplicates(object[["mouseGeneName"]]))
})
