context("transmit : NCBI FTP")

skip_on_appveyor()
skip_on_docker()

## Note that only FTP is currently supported.
remoteDir <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/"

skip_if_not(hasInternet(remoteDir))

test_that("Get README file", {
    object <- transmit(
        remoteDir = remoteDir,
        pattern = "^README\\.txt$",
        compress = FALSE
    )
    expected <- file.path(getwd(), "README.txt")
    names(expected) <- "README.txt"
    expect_identical(object, expected)
    ## Check that function skips on existing files.
    expect_message(
        object = transmit(
            remoteDir = remoteDir,
            pattern = "^README\\.txt$",
            compress = FALSE
        ),
        regexp = "Skipped"
    )
    unlink(object)
})

test_that("Rename and compress", {
    object <- transmit(
        remoteDir = remoteDir,
        pattern = "^README\\.txt$",
        rename = "readme.txt",
        compress = TRUE
    )
    expected <- file.path(getwd(), "readme.txt.gz")
    names(expected) <- "README.txt"
    expect_identical(object, expected)
    unlink("readme.txt.gz")
})

test_that("Invalid parameters", {
    expect_error(
        object = transmit(
            remoteDir = remoteDir,
            pattern = "XXX"
        ),
        regexp = "match"
    )
    expect_error(
        object = transmit(
            remoteDir = remoteDir,
            pattern = "^README$",
            rename = c("XXX", "YYY")
        ),
        regexp = "match"
    )
    ## Currently only FTP is supported.
    expect_error(
        object = transmit(
            remoteDir = "http://steinbaugh.com",
            pattern = "^README$"
        ),
        regexp = "ftp"
    )
    expect_error(
        object = transmit(
            remoteDir = "ftp://ftp.wormbase.org/pub/",
            pattern = "^README$"
        ),
        regexp = "remoteFiles"
    )
})
