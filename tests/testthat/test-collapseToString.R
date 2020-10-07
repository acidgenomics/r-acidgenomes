context("collapseToString")

test_that("atomic", {
    expect_identical(
        object = collapseToString(
            object = groceries,
            sort = TRUE,
            removeNA = FALSE,
            unique = FALSE
        ),
        expected = "eggs, eggs, milk, veggies, NA, NA"
    )
    expect_identical(
        object = collapseToString(
            object = groceries,
            sort = TRUE,
            removeNA = TRUE,
            unique = FALSE
        ),
        expected = "eggs, eggs, milk, veggies"
    )
    expect_identical(
        object = collapseToString(
            object = groceries,
            sort = TRUE,
            removeNA = TRUE,
            unique = TRUE
        ),
        expected = "eggs, milk, veggies"
    )
})

test_that("data.frame", {
    expect_identical(
        object = datasets::mtcars %>%
            head() %>%
            collapseToString(sort = TRUE, unique = TRUE) %>%
            .[, "mpg", drop = TRUE],
        expected = mpgString
    )
})

test_that("DataFrame", {
    expect_identical(
        object = mtcars %>%
            as("DataFrame") %>%
            head() %>%
            collapseToString(sort = TRUE, unique = TRUE) %>%
            .[["mpg"]],
        expected = mpgString
    )
})

test_that("integer", {
    expect_identical(
        object = collapseToString(seq_len(5L)),
        expected = "1, 2, 3, 4, 5"
    )
})

test_that("logical", {
    expect_identical(
        object = collapseToString(c(TRUE, FALSE), sort = TRUE),
        expected = "FALSE, TRUE"
    )
    ## NA and NaN should stay fixed even when sorting is enabled
    expect_identical(
        object = collapseToString(c(NA, NaN), sort = TRUE),
        expected = "NA, NaN"
    )
    expect_identical(
        object = collapseToString(c(NaN, NA), sort = TRUE),
        expected = "NaN, NA"
    )
})

test_that("matrix", {
    expect_identical(
        object = mtcars %>%
            as("matrix") %>%
            head() %>%
            collapseToString(sort = TRUE, unique = TRUE) %>%
            .[1L, "mpg", drop = TRUE] %>%
            as.character(),
        expected = mpgString
    )
})

test_that("numeric", {
    expect_identical(
        object = collapseToString(c(3.141593, 6.0221409e+23)),
        expected = "3.141593, 6.0221409e+23"
    )
})

test_that("scalar early return", {
    expect_identical(collapseToString(1L), 1L)
})
