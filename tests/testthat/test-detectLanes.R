context("detectLanes")

test_that("FASTQ files", {
    files <- paste0(
        "sample1",
        paste0("_R", seq_len(2L)),
        paste0("_L00", seq_len(4L)),
        ".fastq.gz"
    )
    expect_identical(
        object = detectLanes(files),
        expected = seq_len(4L)
    )
})

test_that("Detection failure", {
    expect_identical(
        object = detectLanes(c("aaa", "bbb")),
        expected = integer()
    )
})
