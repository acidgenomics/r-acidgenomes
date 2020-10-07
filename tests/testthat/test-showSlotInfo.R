context("showSlotInfo")

test_that("showSlotInfo", {
    expect_output(
        object = showSlotInfo(list(filtered = TRUE)),
        regexp = "filtered: TRUE"
    )
})
