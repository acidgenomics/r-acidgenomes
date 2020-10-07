context("mutate")

x <- as(mtcars, "DataFrame")

test_that("mutateAll", {
    expect_identical(x[["gear"]][1L], 4)  # nolint
    x <- mutateAll(x, fun = log, base = 2L)
    expect_s4_class(x, "DataFrame")
    ## Check that log2 calculation was applied correctly.
    expect_identical(x[["gear"]][1L], 2)  # nolint
    expect_true(hasRownames(x))
})

## > test_that("mutateAll : List column", {
## >     object <- DataFrame(a = I(list(c(1L, 2L), c(3L, 4L))))
## >     expect_s4_class(object, "DataFrame")
## >     expect_identical(
## >         object = dim(object),
## >         expected = c(2L, 1L)
## >     )
## >     expect_identical(
## >         object = mutateAll(object, fun = toString),
## >         expected = DataFrame(a = c("1, 2", "3, 4"))
## >     )
## > })

test_that("mutateAt", {
    x <- mutateAt(x, vars = c("mpg", "cyl"), log, base = 2L)
    expect_s4_class(x, "DataFrame")
    expect_true(hasRownames(x))
})

## > test_that("mutateAt : List column", {
## >     object <- DataFrame(
## >         a = I(list(c(1L, 2L), c(3L, 4L))),
## >         b = I(list(NULL, NULL))
## >     )
## >     expect_s4_class(object, "DataFrame")
## >     expect_identical(
## >         object = dim(object),
## >         expected = c(2L, 2L)
## >     )
## >     expect_identical(
## >         object = mutateAt(object, vars = "a", fun = toString),
## >         expected = DataFrame(
## >             a = c("1, 2", "3, 4"),
## >             b = I(list(NULL, NULL))
## >         )
## >     )
## > })

test_that("mutateIf", {
    x <- mutateIf(x, predicate = is.double, fun = as.integer)
    expect_s4_class(x, "DataFrame")
    expect_true(hasRownames(x))
})

test_that("transmuteAt", {
    x <- transmuteAt(x, vars = c("mpg", "cyl"), log, base = 2L)
    expect_s4_class(x, "DataFrame")
    expect_true(hasRownames(x))
})

test_that("transmuteIf", {
    x <- transmuteIf(x, predicate = is.double, fun = as.integer)
    expect_s4_class(x, "DataFrame")
    expect_true(hasRownames(x))
})
