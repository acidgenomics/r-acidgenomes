context("minimalSampleData")

test_that("minimalSampleData", {
    expect_identical(
        object = minimalSampleData(c("sample 1", "sample 2")),
        expected = DataFrame(
            sampleName = factor(c("sample 1", "sample 2")),
            row.names = factor(c("sample_1", "sample_2"))
        )
    )
})
