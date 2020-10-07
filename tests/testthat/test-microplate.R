context("microplate")

test_that("96-well plate format", {
    plate <- microplate(plates = 1L, wells = 96L)
    expect_is(plate, "character")
    expect_identical(length(plate), 96L)
    expect_identical(
        head(plate),
        c("1-A01", "1-A02", "1-A03", "1-A04", "1-A05", "1-A06")
    )
    expect_identical(
        tail(plate),
        c("1-H07", "1-H08", "1-H09", "1-H10", "1-H11", "1-H12")
    )
})

test_that("384-well plate format", {
    plate <- microplate(plates = 1L, wells = 384L)
    expect_is(plate, "character")
    expect_identical(
        tail(plate),
        c("1-P19", "1-P20", "1-P21", "1-P22", "1-P23", "1-P24")
    )
})

test_that("Multiple plates", {
    expect_identical(
        object = length(microplate(plates = 2L, wells = 96L)),
        expected = 192L
    )
    expect_identical(
        object = length(microplate(plates = 2L, wells = 384L)),
        expected = 768L
    )
})

test_that("Control wells", {
    plate <- microplate(controls = 3L)
    expect_identical(plate[[1L]], "1-A04")
})

test_that("Prefix", {
    plates <- microplate(prefix = "cherrypick")
    expect_identical(plates[[1L]], "cherrypick-1-A01")
})

test_that("Invalid parameters", {
    expect_error(
        object = microplate(plates = -1L),
        regexp = "isPositive"
    )
    expect_error(
        object = microplate(wells = 4L),
        regexp = "isSubset"
    )
    expect_error(
        object = microplate(controls = -1L),
        regexp = "isInRange"
    )
    expect_error(
        object = microplate(prefix = c("a", "b")),
        regexp = "isString"
    )
})
