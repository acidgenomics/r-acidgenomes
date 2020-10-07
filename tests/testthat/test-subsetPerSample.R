SCE <-  # nolint
    structure("SingleCellExperiment", package = "SingleCellExperiment")



context("subsetPerSample")

test_that("List mode", {
    x <- subsetPerSample(sce, assignAndSave = FALSE)
    expect_is(x, "list")
    expect_identical(
        object = names(x),
        expected = paste0("sample", seq_len(2L))
    )
    expect_identical(
        lapply(x, class),
        list(
            sample1 = SCE,
            sample2 = SCE
        )
    )
})

## This is useful for larger datasets.
test_that("Assign and save mode", {
    env <- new.env()
    files <-
        subsetPerSample(
            object = sce,
            assignAndSave = TRUE,
            envir = env,
            dir = "subsetPerSample"
        )
    expect_identical(
        object = files,
        expected = c(
            sample1 = realpath(file.path("subsetPerSample", "sample1.rds")),
            sample2 = realpath(file.path("subsetPerSample", "sample2.rds"))
        )
    )
    expect_identical(
        object = ls(env),
        expected = paste0("sample", seq_len(2L))
    )
    expect_identical(
        object = sort(list.files("subsetPerSample")),
        expected = paste0("sample", seq_len(2L), ".rds")
    )
    unlink("subsetPerSample", recursive = TRUE)
})
