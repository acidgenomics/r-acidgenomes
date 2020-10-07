context("multiassignAsEnvir")

test_that("multiassignAsEnvir", {
    object <- multiassignAsEnvir(rse, sce, envirName = "example")
    expected <- c("rse", "sce")
    expect_identical(object, expected)
    expect_identical(
        object = ls(example),
        expected = expected
    )
    expect_error(
        object = multiassignAsEnvir(rse, envirName = parent.frame()),
        regexp = "isString.*envirName"
    )
})
