context("removeNA")

## Support for vectors (using `stats::na.omit`).
## This will return structure attributes about original size, with class omit.

test_that("removeNA", {
    mapply(
        object = list(
            character = c("hello", "world", NA),
            numeric = c(1L, 2L, NA),
            DataFrame = DataFrame(
                a = c("A", NA, "C"),
                b = c(NA, NA, NA),
                c = c("B", NA, "D"),
                row.names = c("x", "y", "z")
            )
        ),
        expected = list(
            character = structure(
                .Data = c("hello", "world"),
                na.action = structure(3L, class = "omit")
            ),
            numeric = structure(
                .Data = c(1L, 2L),
                na.action = structure(3L, class = "omit")
            ),
            DataFrame = DataFrame(
                a = c("A", "C"),
                c = c("B", "D"),
                row.names = c("x", "z")
            )
        ),
        FUN = function(object, expected) {
            expect_identical(
                object = removeNA(object),
                expected = expected
            )
        },
        SIMPLIFY = FALSE
    )
})
