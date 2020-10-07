context("detectHPC")

test_that("Local machine", {
    expect_identical(detectHPC(), FALSE)
})

test_that("LSF", {
    Sys.setenv("LSF_ENVDIR" = "XXX")
    expect_identical(detectHPC(), "LSF")
    Sys.unsetenv("LSF_ENVDIR")
})

test_that("SLURM", {
    Sys.setenv("SLURM_CONF" = "XXX")
    expect_identical(detectHPC(), "SLURM")
    Sys.unsetenv("SLURM_CONF")
})
